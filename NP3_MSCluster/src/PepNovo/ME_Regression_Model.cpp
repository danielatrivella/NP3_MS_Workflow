#include "ME_REG.h"

/************************************************************************************
 Takes all the samples in the dataset with the right label and calculates the probs
 sorts them. Let x be the probability at the desired percentile, and t be the
 target probability. The function retutns y, s.t. x^y=t
 This is used as a scaling mechanism to bring the probabilities to a desired form
 without making violoations (values<0 || >1).
*************************************************************************************/
float ME_Regression_Model::calc_log_scaling_constant(int label, 
									const ME_Regression_DataSet& ds, 
									float target_prob) const
{
	vector<float> probs;

	int i;
	for (i=0; i<ds.samples.size(); i++)
	{
		if (ds.samples[i].label != label)
			continue;

		probs.push_back(p_y_given_x(label,ds.samples[i]));
	}

	sort(probs.begin(),probs.end());
	double sum_logs = 0;
	for (i=0; i<probs.size(); i++)
		sum_logs+=log(probs[i]);

	double avg_log = sum_logs/probs.size();
	double log_target = log(target_prob);

	float diff = log_target - avg_log;



/*	int idx = percentile * probs.size();

	float x = probs[idx];
	if (x<=0 || x>=1.0)
		return 1.0;

	float y = log(target_prob)/log(x);*/

	float y=exp(diff);

	if (1)
	{
		cout << "Target prob: " << setprecision(3) << target_prob <<  endl;
		cout << "Mult factor: " <<  y << endl;
		float p;
		for (p=0.05; p<1; p+=0.05)
		{
			float prob = probs[int(p*probs.size())];
			cout << p << "\t" << fixed << setprecision(3) << prob << "\t---->\t" << y*prob << endl;
		}
	}

	return y;
}

// sets the weights to the model returns a constant probability p for the class 0 samples
void ME_Regression_Model::set_weigts_for_const_prob(float p)
{
	int i;
	for (i=0; i<f_weights.size(); i++)
		f_weights[i]=0;

	// set the weight of the constant feature (0)
	float x = log(p);
	f_weights[0]=x/(1.0-x);
}

/************************************************************
Shows what features and values contribute to the exp_sum
*************************************************************/
void ME_Regression_Model::report_exp_sums(const ME_Regression_Sample& sam) const
{
	int i;
	double sum=0;
	for (i=0; i<sam.f_vals.size(); i++)
	{
		int f_idx = sam.f_vals[i].f_idx;
		cout << setw(4) << left << f_idx << " ";
		cout << setw(6) << f_weights[f_idx] << " * " << setw(6) << sam.f_vals[i].val <<
			" = " << setw(8) << f_weights[f_idx] * sam.f_vals[i].val;
		sum+= f_weights[f_idx] * sam.f_vals[i].val;
		cout << "   (" << setw(6) << sum << ")" << endl;
	}
}




void ME_Regression_Model::print_ds_probs(const ME_Regression_DataSet &ds) const
{
    int i;

    vector<float> pr_pos,pr_neg;
    pr_pos.clear();
    pr_neg.clear();
    
    for (i=0; i<ds.samples.size(); i++)
    {
        float prob = p_y_given_x(0,ds.samples[i]);
     
        if (ds.samples[i].label==0)
        {
            pr_pos.push_back(prob);
        }
        else
            pr_neg.push_back(prob);
    }

    sort(pr_pos.begin(),pr_pos.end());
    sort(pr_neg.begin(),pr_neg.end());

    printf("\n\n");
    if (pr_pos.size()<10)
    {
        printf("#pos samples %d.\n",pr_pos.size());
    }
    else
    {
        double av=0;
        for (i=0; i<pr_pos.size(); i++)
            av+=pr_pos[i];
        av /= pr_pos.size();
        printf("#pos samples %d, avg prob = %.3f\n",pr_pos.size(),av);
        
        // prints avgs of tenths of the values
        int ts=pr_pos.size()/10;
        int p=0;
        for (i=0; i<9; i++)
        {
            int next=p+ts;
            int j;
            double av=0;
            for (j=p; j<next; j++)
                av+=pr_pos[j];

            printf("%.4f  ",av/ts);
            p+=ts;
        }

        av=0;
        for (i=p; i<pr_pos.size(); i++)
            av+=pr_pos[i];

        printf("%.4f\n",av/(pr_pos.size()-p));
    }

    if (pr_neg.size()<10)
    {
        printf("#neg samples %d.\n",pr_neg.size());
    }
    else
    {
        double av=0;
        for (i=0; i<pr_neg.size(); i++)
            av+=pr_neg[i];
        av /= pr_neg.size();
        printf("#neg samples %d, avg prob = %.3f\n",pr_neg.size(),av);
        
        // prints avgs of tenths of the values
        int ts=pr_neg.size()/10;
        int p=0;
        for (i=0; i<9; i++)
        {
            int next=p+ts;
            int j;
            double av=0;
            for (j=p; j<next; j++)
                av+=pr_neg[j];

            printf("%.4f  ",av/ts);
            p+=ts;
        }

        av=0;
        for (i=p; i<pr_neg.size(); i++)
            av+=pr_neg[i];

        printf("%.4f\n",av/(pr_neg.size()-p));
    }
}






void ME_Regression_Model::print_ds_histogram(const ME_Regression_DataSet& ds) const
{
    int i;

    vector<float> pr_pos,pr_neg;
    pr_pos.clear();
    pr_neg.clear();
    
    for (i=0; i<ds.samples.size(); i++)
    {
        float prob = p_y_given_x(0,ds.samples[i]);
     
        if (ds.samples[i].label==0)
        {
            pr_pos.push_back(prob);
        }
        else
            pr_neg.push_back(prob);
    }

    sort(pr_pos.begin(),pr_pos.end());
    sort(pr_neg.begin(),pr_neg.end());

	vector<float> counts_p,counts_n;
	counts_p.resize(20);
	counts_n.resize(20);
	for (i=0; i<20; i++)
	{
		counts_p[i]=0;
		counts_n[i]=0;
	}

	for (i=0; i<pr_pos.size(); i++)
		counts_p[(int)(pr_pos[i]*20)]++;

	for (i=0; i<pr_neg.size(); i++)
		counts_n[(int)(pr_neg[i]*20)]++;

	printf("\nRange\t\tPos\tNeg\n");

	float c_pos=0;
	float c_neg=0;
	for (i=0; i<20; i++)
	{
		c_pos+=counts_p[i];
		c_neg+=counts_n[i];
		printf("%.2f - %.2f \t%.3f\t%.3f\n",i*0.05,(i+1)*0.05,c_pos/pr_pos.size(),
			c_neg/pr_neg.size());
	}
    printf("\n\n");   
}


void ME_Regression_Model::write_regression_model(ostream& os) const
{
	int i;
	os << fixed << num_features << endl;
	os << scientific << setprecision(8);
	for (i=0; i<num_features; i++)
		 os << f_weights[i] << endl;
}

void ME_Regression_Model::read_regression_model(istream& is)
{
	char buff[64];
	is.getline(buff,64);

	istringstream iss(buff);
	num_features=-1;
	iss >> num_features;
	if (num_features<0)
	{
		cout << "Error reading ME regression model: " << buff << endl;
		exit(1);
	}

	f_weights.clear();

	if (num_features>10000)
	{
		cout << "Warning, model with too many features: " << num_features << endl;
	}

	if (num_features==0)
	{
		cout << "Warning, 0 features in ME!" << endl;
		has_weights=false;
		return;
	}

	f_weights.resize(num_features,0);

	int i;
	for (i=0; i<num_features; i++)
	{
		is.getline(buff,32);
		istringstream iss(buff);
		iss >> f_weights[i];
	}

	num_classes=2;
	has_weights=true;

	for (i=0; i<f_weights.size(); i++)
		if (f_weights[i] !=0)
			break;
	
	// ignore models where all weights are 0 (these are bad models that have not converged)
	if (i==f_weights.size())
	{
		f_weights.clear();
		has_weights = false;
	}
}

