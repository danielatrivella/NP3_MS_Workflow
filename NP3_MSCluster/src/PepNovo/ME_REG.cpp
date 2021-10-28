#include "ME_REG.h"



/*********************************************************
//
// Vector and matrix multiplication functions.
//
**********************************************************/

inline double vector_self_dot_prod(const vector<double>& v)
{
	double dp=0;
	int i;

	for (i=0; i<v.size(); i++)
	{
		dp += v[i]*v[i];
	}

/*	i=0;
	while (i<v.size())
	{
		const int end_idx = (i + 1000>v.size() ? v.size() : i+1000);
		int j;
		
		double tmp=0;
		for (j=i; j<end_idx; j++)
			tmp += v[i]*v[i];
		dp+=tmp;
		i=end_idx;
	}*/

	return dp;
}


inline double vector_dot_prod(const vector<double>& v1, const vector<double>& v2 )
{
	double dp=0;
	int i;

	if (v1.size() != v2.size())
	{
		printf("mismatch in vector dimensions!\n");
		exit(1);
	}
	

	for (i=0; i<v1.size(); i++)
		dp += v1[i]*v2[i];

/*	i=0;
	while (i<v1.size())
	{
		const int end_idx = (i + 1000>v1.size() ? v1.size() : i+1000);
		int j;
		
		double tmp=0;
		for (j=i; j<end_idx; j++)
			tmp += v1[i]*v2[i];
		dp+=tmp;
		i=end_idx;
	}*/


	return dp;
}

/****************************************************************
Returns the result of multipling a vecotr by matrix. matrix should
come as a vector of rows.
*****************************************************************/
inline void row_times_matrix(const vector<double>& v, 
							 const vector< vector<double> >& m,
							 vector<double>& res)
{
	const int len=m[0].size();
	const int v_len=v.size();
	int i;

	if (v_len != m.size())
	{
		printf("mismatch in vector and matrix dimensions!\n");
		exit(1);
	}
	res.resize(len,0);
	
	for (i=0; i<len; i++)
	{
		double tmp=0;
		int j;
		for (j=0; j<v_len; j++)
			tmp += v[j]*m[j][i];
		res[i]+=tmp;
	}
}




struct fpair {
	fpair(int i, value_t v) : sam_idx(i), val(v) {}
	int sam_idx;
	value_t val;
};



/*********************************************************************/
//
//   Logisitic Regression Functions
//
/*********************************************************************/





/*********************************************************************
Calculates the vector for the first derivatives with the current lambdas
This function assumes that the values of p and fp are already computed.
**********************************************************************/
void ME_Regression_Model::calc_neg_first_derivative(const ME_Regression_DataSet& ds, 
													const vector<double>& p,
													vector<double>& vals) const
{
	int i;

	vals.clear();
	vals.resize(num_features,0);

	for (i=0; i<ds.num_samples; i++)
	{
		int f,label = ds.samples[i].label;
		const ME_Regression_Sample& sam = ds.samples[i];
		const double& sam_prob = p[i];

		if (label == 0)
		{
			double one_minus_p = 1 - sam_prob;
			for (f=0; f<sam.f_vals.size(); f++)
			{
				const int& f_idx = sam.f_vals[f].f_idx;
				vals[f_idx] += sam.weight * ( sam.f_vals[f].val * one_minus_p);
			}
		}
		else
		{
			for (f=0; f<sam.f_vals.size(); f++)
			{
				const int& f_idx = sam.f_vals[f].f_idx;
				vals[f_idx] -= sam.weight * ( sam.f_vals[f].val * sam_prob);
			}
		}
	}
}




/******************************************************************
Calculates the matrix for the Hesian of the current lambdas.
This function assumes that the values of p and fp are already computed.
*******************************************************************/
void ME_Regression_Model::calc_Hesian2(const ME_Regression_DataSet& ds, 
									   const vector<double>& p,
									   vector< vector<double> >& hes) const
{
	int i;
	int num_features=ds.num_features;
	int num_classes=ds.num_classes;

	hes.resize(num_features);
	for (i=0; i<num_features; i++)
	{
		hes[i].clear();
		hes[i].resize(num_features,0);
	}
	

	for (i=0; i<ds.num_samples; i++)
	{
		int f;
		double p_times_one_minus_p = p[i] * ( 1.0 - p[i]);
		const ME_Regression_Sample& sam = ds.samples[i];
		
		for (f=0; f<sam.f_vals.size(); f++)
		{
			const int& f_idx = sam.f_vals[f].f_idx;
			int g;
			for (g=f; g<sam.f_vals.size(); g++)
			{
				const int &g_idx = sam.f_vals[g].f_idx;

				hes[f_idx][g_idx] += sam.weight * sam.f_vals[f].val * 
									 sam.f_vals[g].val * p_times_one_minus_p;
			}
		}
	}

	// fill symmetric values
	int g;
	for (g=0; g<ds.num_features; g++)
	{
		int f;
			
		for (f=g+1 ; f<ds.num_features; f++)
			hes[f][g]=hes[g][f];
	}
}




/************************************************************************
// calculates p(y|x) for all x,y
// first index in p is class index 0..k-1
// second index in p is sample number 0..n-1
*************************************************************************/
void ME_Regression_Model::calc_p_for_all_samples(const ME_Regression_DataSet& ds, 
												 vector<double>& p) const
{
	int i;

	int num_classes = ds.num_classes;
	int num_features = ds.num_features;
	int num_samples = ds.num_samples;

	if (p.size() != num_samples)
	{
		p.resize(num_samples);
	}

	for (i=0; i<num_samples; i++)
		p[i] = this->p_y_given_x(0,ds.samples[i]);
}




/*************************************************************
Calc the log likelihood of the dataset
**************************************************************/
double ME_Regression_Model::log_likelihood(const ME_Regression_DataSet& ds) const
{
	double ll=0;
	int i;

	for (i=0; i<ds.num_samples; i++)
	{
		double p=p_y_given_x(ds.samples[i].label,ds.samples[i]);
		ll+=ds.samples[i].weight*log(p);
	}

//	cout << "Likelihood: " << ll << endl;

	return ll;
}



/***********************************************************************************
Trains the ME model using the cg method.
Returns false if there was no convergence (numerical stability issues)
************************************************************************************/
bool ME_Regression_Model::train_cg(const ME_Regression_DataSet& ds, 
								   int max_iterations,
								   double epsilon,
								   int reset_rounds)
{
	int i,k,f,n;
	double delta_old, delta_new, delta_zero, eps_sqr_delta_zero, prev_likelihood,org_likelihood;
	vector<double> r,d,p;
	
	has_weights=false;
	num_features = ds.num_features;
	num_classes = ds.num_classes;

	i=0;
	k=0;
	n=num_features;

	if (reset_rounds == 0)
		reset_rounds = n-1;

	vector<double> prev_good_weights;
	bool moderate_advance = false;
	prev_good_weights.clear();

RESET_ROUNDS:	
		

	printf("Max iterations: %d  reset CG after %d rounds\n",max_iterations,reset_rounds);

	// initialize weights (use previous values if weights already set)


	if (prev_good_weights.size() == n)
	{
		f_weights = prev_good_weights;
		moderate_advance = true;
	}
	else
	{
		f_weights.clear();
		f_weights.resize(n,0);
	}
	
	
	calc_p_for_all_samples(ds,p);
	calc_neg_first_derivative(ds,p,r);

	d=r;
	delta_new=vector_self_dot_prod(r);
	delta_zero=delta_new;
	eps_sqr_delta_zero = epsilon * epsilon * delta_zero;



	while (i<max_iterations && delta_new>eps_sqr_delta_zero)
	{
		int j=0,j_max = 10;
		double delta_d = vector_self_dot_prod(d);
		double newton_eps_sqr = 0.001;
		double alpha,beta;
		vector< vector<double> > hes;
		vector<double> temp;
		
		prev_good_weights = f_weights;
		
		bool extra_newton_raphson_iter = false;
		do
		{
			double top,denominator;
			int f;
		
			calc_Hesian2(ds,p,hes);

			top = vector_dot_prod(r,d);
			row_times_matrix(d,hes,temp);
			denominator = vector_dot_prod(temp,d);

			alpha = top / denominator; // the minus is in r (which is the negative first derivative)

			if (moderate_advance)
				alpha *= 0.25;

			if (! (alpha == alpha))
			{
				cout << "Alpha NAN!" << endl;
				exit(1);
			}

			// x = x + alpha*d   - just update the weights
			bool halt=false;
			for (f=0; f<num_features; f++)
				f_weights[f] += alpha * d[f];

			// recalculate the first derivative and p and fp
			calc_p_for_all_samples(ds,p);
			calc_neg_first_derivative(ds,p,r);
			j++;

			extra_newton_raphson_iter = (j<j_max && alpha * alpha * delta_d > newton_eps_sqr);
		}
		while(extra_newton_raphson_iter);
		

		delta_old = delta_new;
		delta_new = vector_self_dot_prod(r);
		beta = delta_new / delta_old;

		if (! (beta == beta))
		{
			cout << "Beta NAN!" << endl;
			exit(1);
		}

		for (f=0; f<num_features; f++)
			d[f] = r[f] + beta * d[f];

		k++;

		if (k == reset_rounds || vector_dot_prod(r,d) <= 0)
		{
			d=r;
			k=0;
		}

		// check if we are having a log likelihood problem 
		// can happen because of numerical stability problems
		double likelihood = log_likelihood(ds);
		if (i>0)
		{
			double like_diff = prev_likelihood - likelihood;
			if (like_diff >0)
			{
				cout << i << "\t" << "diff: " << like_diff << endl;
				cout << i << "\t" << "prev: " << prev_likelihood << endl;
				cout << i << "\t" << "curr: " << likelihood << endl;

				if (fabs(like_diff/prev_likelihood)>0.0001)
				{
					reset_rounds = reset_rounds /2;
					if (reset_rounds<1)
					{
						// maybe try using IIS for a few rounds, and use those
						// weights as a starting position for the CG...

						printf("Couldn't train model!\n");
						f_weights.clear();
						f_weights.resize(num_features,0);
						return false;
					}


					if (i>25 || likelihood/org_likelihood<0.5)
					{
						printf("Found likelihood decrease, not training this model further!\n");
						break;
					}

					printf("Likelihood decrease detected, reducing reset rounds to %d\n",reset_rounds);
					i=0;
					goto RESET_ROUNDS;
				}
			}
			else
				moderate_advance=false;

			// means we have a problem
			if (likelihood < -9E200)
			{
				// maybe try using IIS for a few rounds, and use those
				// weights as a starting position for the CG...

				f_weights.clear();
				f_weights.resize(num_features,0);

				printf("Couldn't train model!\n");
				return false;
			}

		}
		else
			org_likelihood = likelihood;

		if (i % 10 == 0)
			printf("%d  NR iterations: %d   Log-likelihood: %.8f\n",i,j,likelihood);
		prev_likelihood = likelihood;
		i++;
	}

	has_weights=true;
	return true;
}

