#include "ME_REG.h"
#include "../Common/auxfun.h"



/*****************************************************************/
//
//
//		LOGISTIC REGRESSION
//
/*****************************************************************/


void ME_Regression_DataSet::add_sample(const ME_Regression_Sample& sam)
{
	if (sam.f_vals.size()>0)
	{
		samples.push_back(sam);
		num_samples=samples.size();
	}
}

void ME_Regression_DataSet::add_samples(const vector<ME_Regression_Sample>& new_samples)
{
	int i;
	for (i=0; i<new_samples.size(); i++)
		add_sample(new_samples[i]);
	
}


// calculate statistics when all samples are present:
// f#, weight of samples in each class
void ME_Regression_DataSet::tally_samples(bool print)
{
	int i;

	num_samples = samples.size();

	// count samples from each class
	class_weights.resize(num_classes);
	for (i=0; i<num_classes; i++)
		class_weights[i]=0;

	
	for (i=0; i<num_samples; i++)
		class_weights[samples[i].label]+=samples[i].weight;		

	
	total_weight=0;
	for (i=0; i<num_classes; i++)
		total_weight+=class_weights[i];

}


// checks that features are within bounds and values are ok (not nan)
int ME_Regression_DataSet::check_samples(bool remove_bad)
{
	int num_bad=0;
	int i;

	for (i=0; i<samples.size(); i++)
	{
		bool is_ok=true;
		vector<fval>& fvals = samples[i].f_vals;
		int j;
		for (j=0; j<fvals.size(); j++)
		{
			if (fvals[j].f_idx<0 || fvals[j].f_idx>=num_features)
				is_ok=false;
			if (! (fvals[j].val == fvals[j].val))
			{
				cout << "Warning! bad value (nan or #INF)" << endl;
				is_ok=false;
			}
		}
		if (! is_ok)
		{
			num_bad++;
			if (remove_bad)
			{
				samples[i]=samples[samples.size()-1];
				samples.pop_back();
			}
		}
	}
	return num_bad;
}




// output dataset, if null goes to screnn
void ME_Regression_DataSet::write_data_to_file(const char *file) const
{
	FILE *stream;
	int i;

	if (file)
	{
		stream=fopen(file,"w");
		if (! stream)
		{
			printf("Couldn't open %s for writing!\n",file);
			exit(1);
		}
	}
	else
		stream=stdout;

	for (i=0; i<num_samples; i++)
	{
		int f;
		fprintf(stream,"CLASS %d %g %d\n",samples[i].label,samples[i].weight,
			static_cast<int>(samples[i].f_vals.size()));

		for (f=0; f<num_features-1; f++)
			fprintf(stream,"%d %g ",samples[i].f_vals[f].f_idx,samples[i].f_vals[f].val);
		
		fprintf(stream,"%d %g\n",samples[i].f_vals[f].f_idx,samples[i].f_vals[f].val);
	}

	if (file)
		fclose(stream);
}


// Splits the dataset into two new mutually exclusive set, must supply
// a vector with the indices of the test set, all the other samples
// are sent to the trainig set.
void ME_Regression_DataSet::split_set(ME_Regression_DataSet& training,ME_Regression_DataSet& test, 
						   vector<int> test_idxs) const
{
	int i;

	vector<int> taken;

	training.clear();
	test.clear();

	training.samples.reserve(test_idxs.size());
	test.samples.reserve(num_samples-test_idxs.size());

	training.num_classes=num_classes;
	training.num_features=num_features;

	test.num_classes=num_classes;
	test.num_features=num_features;

	taken.resize(num_samples);
	for (i=0; i<num_samples; i++)
		taken[i]=0;

	for (i=0; i<test_idxs.size(); i++)
	{
		if (test_idxs[i]>=num_samples || test_idxs[i]<0 )
		{
			printf("Illegal index in data set split: %d\n",test_idxs[i]);
			exit(1);
		}

		taken[test_idxs[i]]=1;
		test.add_sample(samples[test_idxs[i]]);
	}

	for (i=0; i<num_samples; i++)
		if (! taken[i])
			training.add_sample(samples[i]);


	training.tally_samples();
	test.tally_samples();

}


// calibrates class weights so class 0 has the desired relative weight
void ME_Regression_DataSet::calibrate_class_weights(double class_0_weight)
{
	tally_samples();
	if (class_weights[0] == 0 || class_weights[1] == 0)
	{
		printf("Error: original class weights are %g %g\n",class_weights[0], class_weights[1]);
		exit(1);
	}

	int i;
	double mult0=class_0_weight / (class_weights[0]/total_weight);
	double mult1=(1.0 - class_0_weight) / (class_weights[1]/total_weight);

	for (i=0; i<samples.size(); i++)
		if (samples[i].label == 0)
		{
			samples[i].weight *= mult0;
		}
		else
			samples[i].weight *= mult1;

	tally_samples();
}


// calibrates the weight of the entire samples so the new weight adds
// up to *total_weight*
void ME_Regression_DataSet::rescale_dataset_weights(double new_total_weight)
{
	double mult = new_total_weight / total_weight;
	int i;

	for (i=0; i<num_samples; i++)
		samples[i].weight *= mult;

	tally_samples();
}


/*********************************************************************
 changes the weights of the data points so according to the weights
 given in ratios (must have k ratios for values of prop that
 are assumed to be 0,1,..,k-1
**********************************************************************/
void ME_Regression_DataSet::rescale_dataset_according_to_prop(const vector<double>& ratios)
{
	int i;
	double total=0;
	for (i=0; i<ratios.size(); i++)
		total+=ratios[i];

	vector<double> norm_ratios;
	norm_ratios.resize(ratios.size());
	for (i=0; i<ratios.size(); i++)
		norm_ratios[i]=ratios[i]/total;

	// collect info on current ratios
	vector<double> curr_ratios;
	curr_ratios.resize(ratios.size(),0);

	for (i=0; i<num_samples; i++)
		curr_ratios[samples[i].prop] += samples[i].weight;

	total=0;
	for (i=0; i<curr_ratios.size(); i++)
		total+=curr_ratios[i];

	for (i=0; i<curr_ratios.size(); i++)
		curr_ratios[i] /= total;

	vector<double> mult_vals;
	mult_vals.resize(curr_ratios.size(),0);
	for (i=0; i<mult_vals.size(); i++)
		if (curr_ratios[i]>0)
			mult_vals[i] = norm_ratios[i] / curr_ratios[i];

	// update weights
	for (i=0; i<num_samples; i++)
		samples[i].weight *= mult_vals[samples[i].prop];

	tally_samples();
}

// removes features that have low occurence counts (to avoid over-fitting
// and model stability issues)
void ME_Regression_DataSet::purge_low_count_features(int min_num_occurrences)
{
	vector< vector<int> > feature_occurrences;

	feature_occurrences.resize(num_classes);
	int c;
	for (c=0; c<num_classes; c++)
		feature_occurrences[c].resize(num_features,0);

	int i;
	for (i=0; i<samples.size(); i++)
	{
		const int label = samples[i].label;
		int f;
		for (f=0; f<samples[i].f_vals.size(); f++)
		{
			const int f_idx = samples[i].f_vals[f].f_idx;
			feature_occurrences[label][f_idx]++;
		}
	}

	int f;
	for (f=0; f<num_features; f++)
	{
		int c;
		for (c=0; c<num_classes; c++)
			if (feature_occurrences[c][f]>0)
				break;

		if (c==num_classes) // all classes have 0 occurrences
			continue;
		
		for (c=0; c<num_classes; c++)
			if (feature_occurrences[c][f]<min_num_occurrences)
				break;

		if (c==num_classes)
			continue;

		cout << "PURGING " << f << "\t" << feature_occurrences[0][f] << "\t" << feature_occurrences[1][f] << endl;

		int i;
		for (i=0; i<samples.size(); i++)
			samples[i].remove_feature(f);
	}

	tally_samples();
}


void ME_Regression_DataSet::randomly_remove_samples_with_activated_feature(int label, 
											int feature_idx, float prob_remove)
{
	const int org_size = samples.size();
	int num_removed=0;
	int i;
	int num_bad = 0;
	for (i=0; i<samples.size(); i++)
	{
		if (samples[i].label != label)
			continue;

		if (samples[i].get_feature_value(feature_idx) == 0)
			continue;

		num_bad++;
		
		if 	(myRandom() < prob_remove)
		{
			samples[i]=samples[samples.size()-1];
			samples.pop_back();
			num_removed++;
		}
	}
	cout << "Removed " << num_removed << "/" << num_bad << endl;
	cout << "Original size " << org_size << ", now " << samples.size() << endl;

	tally_samples();
}

// calculates for each feature (and class) the proportion of samples (weight)
// for which the feature has a non-zero value
void ME_Regression_DataSet::calc_feature_non_zero_weights(vector< vector<double> >& ratios,
														  vector< vector<double> >& avg_nz) const
{
	int s,i,j;

	ratios.resize(num_features);
	avg_nz.resize(num_features);
	for (i=0; i<num_features; i++)
	{
		ratios[i].resize(num_classes);
		avg_nz[i].resize(num_classes);
		for (j=0; j<num_classes; j++)
		{
			ratios[i][j]=0;
			avg_nz[i][j]=0;
		}
	}

	for (s=0; s<num_samples; s++)
	{
		int i;
		for (i=0; i<samples[s].f_vals.size(); i++)
		{
			const int f_idx = samples[s].f_vals[i].f_idx;
			if (f_idx>=num_features)
			{
				cout << "Error: sample has feature with index " << f_idx <<  " (max allowed is " << 
					num_features-1 << ")" << endl;
				exit(1);
			}

			if (samples[s].weight>0)
			{
				ratios[f_idx][samples[s].label]+= samples[s].weight;
				avg_nz[f_idx][samples[s].label]+= samples[s].weight * samples[s].f_vals[i].val;
			}
		}
	}

	// normalize according to weights
	for (i=0; i<num_features; i++)
		for (j=0; j<num_classes; j++)
			if (avg_nz[i][j] != 0)
				avg_nz[i][j] /= ratios[i][j];


	for (i=0; i<num_features; i++)
		for (j=0; j<num_classes; j++)
			ratios[i][j] /= class_weights[j];
}



// returns the relative weight in the class of a certain feature
double ME_Regression_DataSet::get_relative_weight_of_feature(int label, int feature_idx) const
{
	double class_weight=0;
	double feature_weight=0;

	int i;
	for (i=0; i<samples.size(); i++)
	{
		if (samples[i].label != label)
			continue;

		class_weight += samples[i].weight;
		if (samples[i].get_feature_value(feature_idx) != 0)
			feature_weight+=samples[i].weight;
	}
	return (feature_weight/class_weight);
}


// sets the weights of samples in the class in such a way that the relative weight of samples
// with non-zero values for the given feature is given in the relative_weight
void ME_Regression_DataSet::scale_samples_to_feature_relative_weight(int label, 
										int feature_idx, double relative_weight)
{
	double class_weight=0;
	double feature_weight=0;
	vector<bool> sams_ind;

	sams_ind.resize(samples.size(),false);
	int i;
	for (i=0; i<samples.size(); i++)
	{
		if (samples[i].label != label)
			continue;

		class_weight += samples[i].weight;
		if (samples[i].get_feature_value(feature_idx) != 0)
		{
			feature_weight+=samples[i].weight;
			sams_ind[i]=true;
		}
	}


	double org_weight=feature_weight/class_weight;

	if (org_weight<=0)
		return;

	double mult_feature = relative_weight/org_weight;
	if (mult_feature<0.2)
		mult_feature=0.2;
	if (mult_feature>5)
		mult_feature=5;

	relative_weight = mult_feature*org_weight;

	double mult_others = (1.0-relative_weight)/(1.0-org_weight);

	for (i=0; i<samples.size(); i++)
		samples[i].weight *= (sams_ind[i] ? mult_feature : mult_others);
	
	tally_samples();
}


/**************************************************************************
Tries to scale 
***************************************************************************/
void  ME_Regression_DataSet::serial_scale(const vector<int>& feature_idxs)
{
	int i;
	for (i=0; i<feature_idxs.size(); i++)
	{
		const int f_idx = feature_idxs[i];
		double ratio = get_relative_weight_of_feature(0,f_idx);
		scale_samples_to_feature_relative_weight(1,f_idx,ratio);
	}
}


void ME_Regression_DataSet::report_feature_statistics(int f_idx, const char *name) const
{
	vector<double> vals0,vals1;
	double avg_nz0=0, avg_nz1=0;
	double wnz0=0, wnz1=0, wz0=0, wz1=0;
	int i;

	for (i=0; i<samples.size(); i++)
	{
		double val = 0;
		int j;
		for (j=0; j<samples[i].f_vals.size(); j++)
		{
			if (samples[i].f_vals[j].f_idx == f_idx)
			{
				val = samples[i].f_vals[j].val;
				break;
			}
		}

		double weight = samples[i].weight;
		int label = samples[i].label;

		if (val != 0)
		{
			if (label == 0)
			{
				wnz0+= weight;
				avg_nz0 += weight * val;
				vals0.push_back(val);
			}
			else
			{
				wnz1+= weight;
				avg_nz1 += weight * val;
				vals1.push_back(val);
			}
		}
		else
		{
			if (label == 0)
			{
				wz0+=weight;
			}
			else
				wz1+=weight;
		}
	}

	if (avg_nz0 != 0)
		avg_nz0/=wnz0;
	if (avg_nz1 != 0)
		avg_nz1/=wnz1;
	
	printf("Statistics for feature %d ",f_idx);
	if (name)
		printf(" %s",name);
	printf("\n");
	printf("Class 0:\n");
	printf("weight samples with non-zero vals: %.3f (%.2f)  samples with zero val: %.3f (%.2f)\n",
			wnz0,wnz0/(wnz0+wz0),wz0,wz0/(wnz0+wz0));

	printf("Avg weighted: %g    non-weighted vals:\n",avg_nz0);
	sort(vals0.begin(),vals0.end());
	
	// prints avgs of tenths of the values
	int ts=vals0.size()/10;
	int p=0;
	for (i=0; i<9; i++)
	{
		int next=p+ts;
		int j;
		double av=0;
		for (j=p; j<next; j++)
			av+=vals0[j];

		printf("%.4f  ",av/ts);
		p+=ts;
	}

	double av=0;
	for (i=p; i<vals0.size(); i++)
		av+=vals0[i];

	printf("%.4f\n",av/(vals0.size()-p));


	printf("Class 1:\n");
	printf("weight samples with non-zero vals: %.3f (%.2f)  samples with zero val: %.3f (%.2f)\n",
			wnz1,wnz1/(wnz1+wz1),wz1,wz1/(wnz1+wz1));

	printf("Avg weighted: %g    non-weighted vals:\n",avg_nz1);
	sort(vals1.begin(),vals1.end());
	
	// prints avgs of tenths of the values
	ts=vals1.size()/10;
	p=0;
	for (i=0; i<9; i++)
	{
		int next=p+ts;
		int j;
		double av=0;
		for (j=p; j<next; j++)
			av+=vals1[j];

		printf("%.4f  ",av/ts);
		p+=ts;
	}

	av=0;
	for (i=p; i<vals1.size(); i++)
		av+=vals1[i];

	printf("%.4f\n\n\n",av/(vals1.size()-p));

}






// extracts all the samples of the given class and puts them in a new dataset
void ME_Regression_DataSet::extract_class_samples(int label, ME_Regression_DataSet& extract) const
{
	int i;

	extract.samples.clear();
	extract.num_samples = 0;
	extract.num_classes = num_classes;

	for (i=0; i<num_samples; i++)
		if (samples[i].label == label)
			extract.add_sample(samples[i]);

	extract.tally_samples();
}


// exctract samples that have a non-zero value for the given feature
void ME_Regression_DataSet::extract_samples_with_activated_feature(int feature_idx,
												ME_Regression_DataSet& extract) const
{
	int i;

	extract.samples.clear();
	extract.num_samples = 0;
	extract.num_classes = num_classes;	

	for (i=0; i<num_samples; i++)
	{
		int j;
		for (j=0; j<samples[i].f_vals.size(); j++)
			if (samples[i].f_vals[j].f_idx == feature_idx && samples[i].f_vals[j].val != 0)
				extract.add_sample(samples[i]);
	}

	extract.tally_samples();
}


// adds the samples from the other dataset, and adjust weights
void ME_Regression_DataSet::add_other_dataset_samples(const ME_Regression_DataSet& other)
{
	int i;

	for (i=0; i<other.num_samples; i++)
		add_sample(other.samples[i]);

	tally_samples();
}





// return all samples in the datatset that have a desired label
void ME_Regression_DataSet::get_samples_with_label(int label, vector<int>& idxs) const
{
	int i;

	idxs.clear();

	for (i=0; i<samples.size(); i++)
		if (samples[i].label== label)
			idxs.push_back(i);
}





// prints info on features (num non zero and p~(f) )
void ME_Regression_DataSet::print_feature_summary(ostream& os, const char **feature_names) const
{
	int i;

	vector< vector<double> > ratios, avg_nz;

	calc_feature_non_zero_weights(ratios,avg_nz);

	for (i=0; i<num_features; i++)
	{
		os << setw(4) << left << i << " ";
		os << setw(10) << left << setprecision(3) << ratios[i][0] << " " << setw(10) << setprecision(3) << left << ratios[i][1] << " ";
		os << " ( " << setw(6) << setprecision(3) << left << avg_nz[i][0] << " , " <<  setw(6) << setprecision(3) << left << avg_nz[i][1] << ") ";
		os << "   " << setw(6);
		if (feature_names)
			cout << feature_names[i];
		cout << endl;
	}
}

void ME_Regression_DataSet::clear(int num_classes)
{

	num_classes=num_classes;  // number of classes k in the data = max_label+1
	num_samples=0;
	num_features=0;
	class_weights.clear();
	if (num_classes>0)
		class_weights.resize(num_classes,0);

	samples.clear();
}



void ME_Regression_Sample::print(const char **feature_names) const
{
	int j;

	if (! feature_names)
	{
		cout << "> " << label << " " <<weight << endl;
		for (j=0; j<f_vals.size(); j++)
			cout << f_vals[j].f_idx << " " << f_vals[j].val << " ";

		cout << endl;
		return;
	}

	cout << "LABEL " << label << ",  weight " << weight << endl;
	for (j=0; j<f_vals.size(); j++)
	{
		cout << f_vals[j].f_idx << "\t" << setprecision(3) << fixed << f_vals[j].val << "\t" <<
			feature_names[f_vals[j].f_idx] << endl;
	}
	cout << endl;

}


void ME_Regression_Sample::remove_feature(int f_idx)
{
	int f;
	for (f=0; f<f_vals.size(); f++)
		if (f_vals[f].f_idx == f_idx)
			break;
	
	if (f==f_vals.size())
		return;

	if (f == f_vals.size()-1)
	{
		f_vals.pop_back();
		return;
	}

	int i;
	for (i=f+1; i<f_vals.size(); i++)
		f_vals[i-1]=f_vals[i];

	f_vals.pop_back();
}

void ME_Regression_DataSet::print() const
{
	int i;

	for (i=0; i<num_samples; i++)
		samples[i].print();
}


void ME_Regression_DataSet::print_summary() const
{
	int j;

	printf("Classes %d\n",num_classes);
	printf("Samples %d\n",num_samples);
	printf("Total weight %.3f\n",total_weight);
	printf("Relative class weights:\n");
	for (j=0; j<num_classes; j++)
		printf("%d - %.4f\n",j,class_weights[j]/total_weight);

}


