//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// histogram.h
//
//   histogram for statistical data
//

#ifndef __HISTOGRAM_H
#define __HISTOGRAM_H

class Histogram {
public:
	Histogram(unsigned _buckets, unsigned _stride, unsigned _minbucket_value)
		: numbuckets(_buckets), stride(_stride), minbucket_value(_minbucket_value)
	{
		maxbucket_value = minbucket_value + ((numbuckets-1)*stride);
		buckets = static_cast<uint64_t *>( calloc(numbuckets, sizeof(uint64_t)) );
		insertions = 0;
		insertions_sum = 0;
		max_value = minbucket_value;
	}

	~Histogram() { free(buckets); }

	void insert(unsigned value)
	{
		++ insertions;
		insertions_sum += value;
		max_value = max(max_value, value);

		if (value < minbucket_value) value = minbucket_value;
		if (value > maxbucket_value) value = maxbucket_value;

		unsigned index = (value - minbucket_value) / stride;
		assert( index < numbuckets );
		++ buckets[index];
	}

	void print(const char *prefix, uint64_t potential_insertions, uint64_t potential_insertions_sum)
	{
		if (prefix != NULL) {
			printf("%s", prefix);
		}

		printf("%3.6g%% :: ", 100.0 * ((double)insertions) / ((double)potential_insertions));
		printf("%3.6g%% :: ", 100.0 * ((double)insertions_sum) / ((double)potential_insertions_sum));

		if (insertions > 0) {
			printf("%3.6g :: ", ((double)insertions_sum / (double)insertions));
			printf("%3.6g :: ", (double)max_value);

			if (((double)buckets[numbuckets-1] / (double)insertions) < 0.2) {
				double sumavg = 0.0;
				for (unsigned i=0; i < numbuckets; ++i) {
					double thisavg = ((i*stride)+minbucket_value) * ((double)buckets[i] / (double)insertions);
					sumavg += thisavg;
				}
				printf("%3.6g ::", sumavg);
			} else {
				printf("NaN ::");
			}

			printf("\n");
			double sumpercent = 0.0;
			for (unsigned i=0; i < numbuckets-1; ++i) {
				double thispercent = 100.0 * ((double)buckets[i] / (double)insertions);
				printf("  [%02u] %3.6g\t%" PRIu64 "", minbucket_value+(i*stride), thispercent, buckets[i]);
				sumpercent += thispercent;
				printf("\n");
			}
			printf("  [%02u] %3.6g", maxbucket_value, 100.0-sumpercent);
		} else {
            printf("0 :: ");
            printf("0 :: ");
            printf("0 :: ");
			printf("\n");
        }
		printf("\n");
	}

private:
	unsigned numbuckets;
	unsigned stride;
	unsigned minbucket_value;
	unsigned maxbucket_value;
	uint64_t *buckets;

	uint64_t insertions;
	uint64_t insertions_sum;
	unsigned max_value;
};


#endif // __HISTOGRAM_H
