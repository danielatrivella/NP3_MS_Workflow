#ifndef __PEAKWEIGHTTABLE_H__
#define __PEAKWEIGHTTABLE_H__

#include "../Common/auxfun.h"

class PeakWeightTable {
public:
	void initWeights(int maxN, float p)
	{
		assert(p>0.0 && p<1.0);
		maxN_ = maxN;
		kWithValueOne_.resize(maxN_+1);
		weights_.resize(maxN_+1);
		for (int n=1; n<=maxN_; n++)
		{
			weights_[n].push_back(0.0);
			vector<float> cdf;
			computeBinomailCDFs(n, p, cdf);
			int i;
			for (i=0; i<=n; i++)
			{
				float w=cdf[i];
				if (w > 0.99)
					break;
				weights_[n].push_back(w);
			}
			kWithValueOne_[n]=i;
		}
	}

	float getWeight(int k, int n)
	{
		while (n>maxN_)
		{
			n >>= 1;
			k >>= 1;
		}

		if (k>=kWithValueOne_[n])
			return 1.0;
		return (weights_[n][k]);
	}

private:
	int maxN_;
	vector<int> kWithValueOne_;
	vector< vector<float> > weights_; // the weights for having k of n peaks
};




#endif

