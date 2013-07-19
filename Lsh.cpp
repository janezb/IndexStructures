#include "StdAfx.h"
#include "Lsh.h"

namespace LshTest
{

// A very simple test where the dat's are sets of integers (TIntV), using the Jaccard distance.
int LshTest()
{
	// keyDesc contains the descriptions of the int-sets in this test.  Each int-set is represented
	// by a sequence of characters in the a..z range; the sets are separated from each other by a space.
	TStr keyDesc = "abcd abce bcef abcde cdghi cd cdefgh";
	//
	typedef TInt TItem;
	typedef TVec<TItem> TDat;
	TStrV keys; keyDesc.SplitOnWs(keys);
	TVec<TDat> dats;
	TLshCollection<TStr, TDat, TMinHasher<TItem>, TJaccardDistance<TItem> > lsh(TMinHasher<TItem>(3, 3));
	for (int i = 0; i < keys.Len(); i++)
	{
		TStr key = keys[i];
		TIntV dat; for (int j = 0; j < key.Len(); j++) dat.Add(key[j] - 'a');
		dat.Sort();
		lsh.Add(key, dat); 
		lsh.Validate();
		dats.Add(dat);
	}
	for (int i = 0; i < keys.Len(); i++)
	{
		TStrV results;
		lsh.ApproxNnQuery(dats[i], 3, results);
		printf("Nn(%s) = [", keys[i].CStr());
		for (int j = 0; j < results.Len(); j++) 
			printf(" %s:%.3f", results[j].CStr(),
				lsh.DistFunc(lsh.GetDat(results[j]), dats[i]));
		printf(" ]\n");
	}
	lsh.Validate();
	return 0;
}


// A class that generates random sets of integers (TIntV's), by picking
// items at random from a universal set that was generated during the call to Init().
class TIntSetGenerator
{
public:
	typedef TInt TItem;
	typedef TVec<TItem> TDat;

	int avgDatSize, universalSetSize;
	TDat universalSet;

	TIntSetGenerator(int UniversalSetSize, int AvgDatSize) : universalSetSize(UniversalSetSize), avgDatSize(AvgDatSize) { }

	void Init(TRnd& rnd)
	{
		universalSet.Clr();
		for (int i = 0, x = 0; i < universalSetSize; i++) {
			x += rnd.GetUniDevInt(1000) + 1;
			universalSet.Add(x); }
		universalSet.Shuffle(rnd);
	}

	void Gen(TRnd &rnd, TDat& dest) const
	{
		const double p = avgDatSize / double(universalSet.Len());
		dest.Clr();
		for (int i = 0; i < universalSet.Len(); i++)
			if (rnd.GetUniDev() < p) dest.Add(universalSet[i]);
		dest.Sort();
	}
};

// Helper functions for use with TSparseVectorGenerator.  GenDimIdx is used to convert 
// a dimension index from 'int' to 'TDimIdx', where the latter may be TInt or TStr.
void GenDimIdx(int dimIdx, TInt& dest) { dest = dimIdx; }
void GenDimIdx(int dimIdx, TStr& dest) { char buf[50]; sprintf(buf, "dim-%d", dimIdx); dest = buf; }

// Helper functions for use with TSparseVectorGenerator.  Convert is used to convert
// a vector from TVec<TKeyDat<TInt, TVal> > into one of the various sparse vector representations:
// TVec<TKeyDat<TDimIdx, TVal> >, TVec<TPair<TDimIdx, TVal> >, THash<TDimIdx, TVal>,
// or even the non-sparse TVec<TVal>.
template<typename TDimIdx, typename TVal> 
void Convert(const TVec<TKeyDat<TInt, TVal> > &src, const TVec<TDimIdx> &dimIndices, TVec<TKeyDat<TDimIdx, TVal> > &dest) { 
	dest.Clr(); for (int i = 0; i < src.Len(); i++) dest.Add(TKeyDat<TDimIdx, TVal>(dimIndices[src[i].Key], src[i].Dat)); }
template<typename TDimIdx, typename TVal> 
void Convert(const TVec<TKeyDat<TInt, TVal> > &src, const TVec<TDimIdx> &dimIndices, THash<TDimIdx, TVal> &dest) { 
	dest.Clr(); for (int i = 0; i < src.Len(); i++) {
		IAssert(! dest.IsKey(dimIndices[src[i].Key]));
		dest.AddDat(dimIndices[src[i].Key], src[i].Dat); }}
template<typename TDimIdx, typename TVal> 
void Convert(const TVec<TKeyDat<TInt, TVal> > &src, const TVec<TDimIdx> &dimIndices, TVec<TVal> &dest) { 
	dest.Clr(); if (src.Empty()) return;
	dest.Gen(src.Last().Key + 1);
	dest.PutAll(TVal(0));
	for (int i = 0; i < src.Len(); i++) dest[src[i].Key] = src[i].Dat; }
template<typename TDimIdx, typename TVal> 
void Convert(const TVec<TKeyDat<TInt, TVal> > &src, const TVec<TDimIdx> &dimIndices, TVec<TPair<TDimIdx, TVal> > &dest) { 
	dest.Clr(); for (int i = 0; i < src.Len(); i++) dest.Add(TPair<TDimIdx, TVal>(dimIndices[src[i].Key], src[i].Dat)); }

// This class generates random sparse vectors.  These will be nDims-dimensional vectors
// with an average of nAvgNonz nonzero components.  TDimIdx can be either TInt or TStr (if you
// want string dimension names of the form "dim-123"), and TDat_ can be one of various sparse vector
// representations: TIntFltPrV, TIntFltKdV, TIntFltH, TStrFltPrV, TStrFltKdV, TStrFltH, 
// or even the non-sparse TFltV.
//
// The sparse vectors are generated in the following way: for each dimension, we choose to keep it 0 with
// probability 1 - nAvgNonz/nDims; otherwise, we choose its value uniformly from the range [-1, +1].
// Thus, the resulting vectors are considerably different than e.g. TF-IDF vectors from a document corpus;
// if you want to simulate those, see TSparseDocGenerator below.
template<typename TDimIdx, typename TDat_>
class TSparseVectorGenerator
{
public:
	typedef TVec<TDimIdx> TDimIdxV;
	typedef TFlt TVal;
	typedef TKeyDat<TDimIdx, TVal> TKd;
	typedef TVec<TKd> TKdV;
	typedef TDat_ TDat;
	TVec<TDimIdx> dimIndices;
	int nDims, nAvgNonz;

	TSparseVectorGenerator(int nDims_, int nAvgNonz_) : nDims(nDims_), nAvgNonz(nAvgNonz_) {
		for (int dimIdx = 0; dimIdx < nDims; dimIdx++) { dimIndices.Add(); GenDimIdx(dimIdx, dimIndices.Last()); } }

protected:

	void Gen_(TRnd& rnd, TVec<TKeyDat<TInt, TVal> >& dest) const
	{
		const double p = nAvgNonz / double(nDims);
		dest.Clr();
		for (int i = 0; i < nDims; i++) 
			if (rnd.GetUniDev() < p) dest.Add(TKeyDat<TInt, TVal>(i, rnd.GetUniDev() * 2 - 1));
	}

public:

	void Init(TRnd& rnd) { }

	void Gen(TRnd& rnd, TDat& dat) const {
		TVec<TKeyDat<TInt, TVal> > kd; Gen_(rnd, kd); Convert(kd, dimIndices, dat); }

};

// A helper class that can hold several representations of the same sparse vector.
template<typename TDimIdx_, typename TVal_>
class TMultiSparseVector
{
public:
	typedef TDimIdx_ TDimIdx;
	typedef TVal_ TVal;
	typedef TKeyDat<TDimIdx, TVal> TKd;
	typedef TPair<TDimIdx, TVal> TPr;
	typedef TVec<TKd> TKdV;
	typedef TVec<TPr> TPrV;
	typedef THash<TDimIdx, TVal> TKdH;
	typedef TVec<TVal> TValV;
	typedef TSparseVectorGenerator<TDimIdx, TKdV> TGen;

	TKdV kdv;
	TKdH kdh;
	TPrV prv;
	TValV valv;

	void Gen(TRnd& rnd, TGen& gen) {
		gen.Gen(rnd, kdv);
		Convert(kdv, gen.dimIndices, kdh);
		Convert(kdv, gen.dimIndices, prv);
		Convert(kdv, gen.dimIndices, valv); }
};

// This function creates 'nVectors' random sparse vectors and performs 'nOps' dot products
// between randomly chosen pairs of vectors.  It tests whether the results are the same
// regardless of which representation of sparse vectors is used (vector of pairs, vectors of KeyDat's,
// hash table, or even a dense vector).
template<class TDimIdx, class TVal>
void TestDots(int nDims, int nAvgNonz, int nVectors, int nOps)
{
	typedef TMultiSparseVector<TDimIdx, TVal> TMsv;
	typedef typename TMsv::TGen TGen;
	typedef TVectorDistanceBase<TDimIdx, TVal> TDist;
	TGen gen(nDims, nAvgNonz);
	TVec<TMsv> vecs;
	TRnd rnd(1234);
	for (int i = 0; i < nVectors; i++) {
		vecs.Add();
		vecs.Last().Gen(rnd, gen); }
	double MaxError = 0;
	for (int i = 0; i < nOps; i++)
	{
		int i1 = rnd.GetUniDevInt(nVectors), i2 = rnd.GetUniDevInt(nVectors);
		double xx, yy, xy; TDist::Dots(vecs[i1].kdv, vecs[i2].kdv, xx, xy, yy);
		double xx1, yy1, xy1; TDist::Dots(vecs[i1].kdh, vecs[i2].kdh, xx1, xy1, yy1);
		MaxError = 0;
		MaxError = TFlt::GetMx(MaxError, fabs(xx - xx1), fabs(yy - yy1), fabs(xy - xy1));
		TDist::Dots(vecs[i1].prv, vecs[i2].prv, xx1, xy1, yy1);
		MaxError = TFlt::GetMx(MaxError, fabs(xx - xx1), fabs(yy - yy1), fabs(xy - xy1));
		TDist::Dots(vecs[i1].valv, vecs[i2].valv, xx1, xy1, yy1);
		MaxError = TFlt::GetMx(MaxError, fabs(xx - xx1), fabs(yy - yy1), fabs(xy - xy1));
		if (MaxError > 1e-8) 
			printf("! %g\n", MaxError);
	}
}

// Creates a TLshCollection and performs a mixture of additions, deletions, modifications
// and nearest neighbor queries.
template<class TDatGenerator, class THasher, class TDistFunc>
int LshTest2_Helper(const int nKeys, TDatGenerator &datGenerator, THasher &hasher, TDistFunc &distFunc, bool testSaving)
{
	typedef TStr TKey;
	typedef typename TDatGenerator::TDat TDat;
	typedef TKeyDat<TKey, TDat> TKd;
	typedef TVec<TKd> TKdV;
	typedef TVec<TKey> TKeyV;
	typedef TDistFunc::TDistance TDistance;
	typedef TKeyDat<TDistance, TInt> TDistIntKd;
	typedef TVec<TDistIntKd> TDistIntKdV;

	int checksum = 0;
	printf("\nLshTest2_Helper(%d keys)\n- Generator: %s\n- Hasher: %s\n- Dist func: %s\n",
		nKeys, typeid(TDatGenerator).name(),
		typeid(THasher).name(),
		typeid(TDistFunc).name());

	TRnd rnd(123); 
	datGenerator.Init(rnd);
	
	TKdV keyDats; // the list of all 'nKeys' (key, dat) pairs used in this experiment
	TIntV keysInHash, keysOutside; // lists of indices into 'keyDats', indicating which keys are currently in the LshCollection and which ones aren't
	for (int i = 0; i < nKeys; i++) {
		char buf[50]; sprintf(buf, "key-%d", i); TStr key = buf;
		TDat dat; datGenerator.Gen(rnd, dat);
		keyDats.Add(TKd(key, dat)); 
		keysOutside.Add(i); }

	//TLshCollection<TStr, TDat, TMinHasher<TItem>, TJaccardDistance<TItem> > lsh(TMinHasher<TItem>(nBands, nFuncPerBand));
	typedef TLshCollection<TKey, TDat, THasher, TDistFunc> TLsh;
	TLsh lsh(hasher, distFunc);

	int nPhases = 10, overallQueryNo = 0;
	for (int phaseNo = 0; phaseNo < nPhases; phaseNo++)
	{
		// Choose a target collection size for this phase; we'll still do both additions and deletions,
		// but the operation that moves the collection towards the target size rather than away from it
		// will be twice as likely.
		int goalSize;
		if (phaseNo == 0) goalSize = rnd.GetUniDevInt((9 * nKeys + 9) / 10, nKeys);
		else if (phaseNo == nPhases) goalSize = 0;
		else goalSize = rnd.GetUniDevInt(nKeys / 10, nKeys - 1);

		while (lsh.Len() != goalSize)
		{
			bool modify = (keysInHash.Len() > 0) && (rnd.GetUniDevInt(3) == 1);
			if (modify)
			{
				// Pick a random key that is currently in the LshCollection and modify its corresponding 'dat'.
				if (false) printf("m");
				int i = rnd.GetUniDevInt(keysInHash.Len());
				int keyIdx = keysInHash[i];
				TKd &kd = keyDats[keyIdx];
				IAssert(lsh.IsKey(kd.Key));
				IAssert(lsh.GetDat(kd.Key) == kd.Dat);
				datGenerator.Gen(rnd, kd.Dat);
				lsh.Add(kd.Key, kd.Dat);
				IAssert(lsh.IsKey(kd.Key));
				IAssert(lsh.GetDat(kd.Key) == kd.Dat);
			}
			else
			{
				int add = rnd.GetUniDevInt(3);
				if (add == 2) add = (lsh.Len() < goalSize ? 1 : 0);
				if (keysOutside.Empty()) add = 0;
				if (keysInHash.Empty()) add = 1;
				if (add)
				{
					if (false) printf("+");
					int i = rnd.GetUniDevInt(keysOutside.Len());
					int keyIdx = keysOutside[i]; keysOutside[i] = keysOutside.Last(); keysOutside.DelLast();
					keysInHash.Add(keyIdx);
					TKd kd = keyDats[keyIdx];
					IAssert(! lsh.IsKey(kd.Key()));
					lsh.Add(kd.Key, kd.Dat);
					IAssert(lsh.IsKey(kd.Key()));
				}
				else // delete
				{
					if (false) printf("-");
					int i = rnd.GetUniDevInt(keysInHash.Len());
					int keyIdx = keysInHash[i]; keysInHash[i] = keysInHash.Last(); keysInHash.DelLast();
					keysOutside.Add(keyIdx);
					TKd kd = keyDats[keyIdx];
					IAssert(lsh.IsKey(kd.Key));
					IAssert(lsh.GetDat(kd.Key) == kd.Dat);
					bool deleted = lsh.Del(kd.Key);
					IAssert(deleted); IAssert(! lsh.IsKey(kd.Key));
				}
			}
			lsh.Validate();
			IAssert(lsh.Len() == keysInHash.Len());
			// Do a nearest-neighbor query.
			for (int queryNo = 0; queryNo < 1; queryNo++, overallQueryNo++)
			{
				if (false && overallQueryNo == 50)
					printf("!");
				int iKey = rnd.GetUniDevInt(nKeys);
				bool qInHash = (iKey < keysInHash.Len());
				int qKeyIdx = (qInHash) ? keysInHash[iKey] : keysOutside[iKey - keysInHash.Len()];
				int qKeyId = lsh.GetKeyId(keyDats[qKeyIdx].Key);
				TKeyV results, results2; const int MaxResults = 10;
				TLsh::TNnQueryReport report;
				int nResults = lsh.ApproxNnQuery(keyDats[qKeyIdx].Dat, MaxResults, 0, &results, &report);
				int nResults2 = lsh.ApproxNnQuery(keyDats[qKeyIdx].Dat, MaxResults, 0, &results2, 0);
				IAssert(nResults == nResults2); 
				for (int i = 0; i < results.Len(); i++) 
					IAssert(results[i] == results2[i]);
				if (qInHash) IAssert(nResults > 0); 
				else if (keysInHash.Empty()) IAssert(nResults > 0);
				else IAssert(nResults >= 0);
				IAssert(nResults == results.Len());
				IAssert(nResults <= MaxResults);
				IAssert(nResults <= keysInHash.Len());
				TIntV resultIdxs;
				// The results should be keys that are currently in the hash table.
				// Also find their indices in 'keyDats'.
				for (int iResult = 0; iResult < nResults; iResult++) {
					int i = 0; while (i < keysInHash.Len() && results[iResult] != keyDats[keysInHash[i]].Key) i++;
					// All the results must of course be keys that are currently in the hash table.
					IAssert(i < keysInHash.Len()); resultIdxs.Add(keysInHash[i]); 
					checksum = (checksum * 13 + i) % 12582917; 
					if (false && overallQueryNo < 200) printf(" %d:%d:%c:%d", overallQueryNo, iResult, 'a' + (checksum % 23), i); }
				// Make sure that the results are all different.
				{TIntV temp = resultIdxs; temp.Sort(); temp.Merge(); IAssert(temp.Len() == resultIdxs.Len());}
				// Make sure the distances are in increasing order.
				if  (qInHash) {
					// The query is in the hash table at the moment, so it should have been
					// discovered as its own nearest neighbor.  So it should be either the first
					// in the list or at least at distance 0 from the first in the list.
					if (resultIdxs[0] != qKeyIdx) {
						TDistance dist = lsh.DistFunc(keyDats[qKeyIdx].Dat, keyDats[resultIdxs[0]].Dat);
						IAssert(dist == 0);
						int i = 0; while (i < resultIdxs.Len() && resultIdxs[i] != qKeyIdx) i++;
						IAssert(i < resultIdxs.Len()); // the query key must have appeared somewhere in the result list 
					}}
				// Ensure that the results were in ascending order.
				for (int i = 1; i < resultIdxs.Len(); i++) {
					TDistFunc::TDistance d1 = lsh.DistFunc(keyDats[qKeyIdx].Dat, keyDats[resultIdxs[i - 1]].Dat),
						d2 = lsh.DistFunc(keyDats[qKeyIdx].Dat, keyDats[resultIdxs[i - 1]].Dat);
					IAssert(d1 <= d2); }
				// Compute the correct results of this query.
				TDistIntKdV correctResults;
				for (int i = 0; i < keysInHash.Len(); i++)
					correctResults.Add(TDistIntKd(lsh.DistFunc(keyDats[qKeyIdx].Dat, keyDats[keysInHash[i]].Dat), keysInHash[i]));
				correctResults.Sort(); 
				TFltV rank; rank.Gen(nKeys); rank.PutAll(-1);
				for (int i = 0; i < correctResults.Len(); )
				{
					int j = i + 1; while (j < correctResults.Len() && correctResults[j].Key == correctResults[i].Key) j++;
					double avgRank = 1 + (i + (j - 1)) / 2.0;
					while (i < j) rank[correctResults[i++].Dat] = avgRank;
				}
				double sumRanks = 0, idealSum = 0;
				int nInter = 0; // the intersection between our results and the correct results (if the former were also cut off at nResults)
				for (int i = 0; i < resultIdxs.Len(); i++) {
					double r = rank[resultIdxs[i]]; IAssert(r > 0);
					if (r <= nResults) nInter++;
					sumRanks += r; idealSum += (i + 1); }
				if (false) {
					printf("qInHash = %s; nResults = %d/%d; sumRanks = %.2f, ideal = %.2f, ratio (bigger is better) = %.4f; intersection = %d, Jaccard score = %.3f\n",
						qInHash ? "yes" : " no", nResults, lsh.Len(), sumRanks, idealSum, idealSum / sumRanks,
						nInter, nInter / float(nResults + nResults - nInter));
					printf("  "); report.Print(); printf("\n"); }
			}
		}
		if (false) lsh.Report();
		if (testSaving)
		{
			PSOut SOut = TMOut::New();
			lsh.Save(*SOut);
			lsh.~TLshCollection();
			memset(&lsh, 0, sizeof(lsh));
			PSIn SIn = ((TMOut *) SOut())->GetSIn();
			new (&lsh) TLsh(*SIn);
		}
	}
	printf("-> checksum = %d\n", checksum);
	return checksum;
}

template<class TDatGenerator, class THasher, class TDistFunc>
void LshTest2_Helper(const int nKeys, TDatGenerator &datGenerator, THasher &hasher, TDistFunc &distFunc)
{
	TDatGenerator dg1 = datGenerator, dg1a = datGenerator, dg2 = datGenerator;
	THasher hasher1 = hasher, hasher1a = hasher, hasher2 = hasher; 
	int chk1 = LshTest2_Helper(nKeys, dg1, hasher1, distFunc, false);
	int chk1a = LshTest2_Helper(nKeys, dg1a, hasher1a, distFunc, false);
	int chk2 = LshTest2_Helper(nKeys, dg2, hasher2, distFunc, true);
	IAssert(chk1 == chk2);
	IAssert(chk1 == chk1a);
}

// This function calls LshTets2_Helper several times to performs several experiments using the same type of dat's
// (namely those  generated by TGen): some of these experiments use the sign finalizer and cosine distance,
// some use the floor finalizer (with the given projectionGranularity) and the Euclidean distance.
template <typename TDimIdx, typename TVal, typename TRandDirStore, typename TGen>
void LshTest2_Helper2(int nKeys, int nBands, int nFuncPerBand, int nDims, int nAvgNonz, double projectionGranularity)
{
	typedef typename TGen::TDat TDat;
	if (true)
	{
		// Sign finalizer, cosine distance
		typedef TRandProjFinalizer_Sign<TVal, TInt> TFinalizer;
		typedef TCosineDistance<TDimIdx, TVal> TDistFunc;
		if (true && nBands <= 8)
			LshTest2_Helper(nKeys, 
				TGen(nDims, nAvgNonz), 
				TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Set<TB8Set>, TRandDirStore, TFinalizer>(nBands, nFuncPerBand, 0, TFinalizer()),
				TDistFunc());
		if (true && nBands <= 32)
			LshTest2_Helper(nKeys, 
				TGen(nDims, nAvgNonz), 
				TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Set<TB32Set>, TRandDirStore, TFinalizer>(nBands, nFuncPerBand, 0, TFinalizer()),
				TDistFunc());
		if (true)
			LshTest2_Helper(nKeys, 
				TGen(nDims, nAvgNonz), 
				TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Set<TBSet>, TRandDirStore, TFinalizer>(nBands, nFuncPerBand, 0, TFinalizer()),
				TDistFunc());
	}
	if (true)
	{
		// Floor finalizer, Euclidean distance
		typedef TInt THashCode;
		typedef TRandProjFinalizer_Floor<TVal, THashCode> TFinalizer;
		typedef TEuclideanDistance<TDimIdx, TVal> TDistFunc;
		if (true)
			LshTest2_Helper(nKeys, 
				TGen(nDims, nAvgNonz), 
				TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Vec<THashCode>, TRandDirStore, TFinalizer>(nBands, nFuncPerBand, nDims, TFinalizer(projectionGranularity)),
				TDistFunc());
	}
}

// This function calls LshTest2_Helper2 several times to perform several experiments with
// several different representations of basically the same kind of sparse vectors.  In all
// cases, the vectors consist of (TDimIdx, TVal) pairs, and in practice these can be represented
// either by a vector of TKeyDat's, a vector of TPair's, or a hash table.
template <typename TDimIdx, typename TVal, typename TRandDirStore>
void LshTest2_Helper3(int nKeys, int nBands, int nFuncPerBand, int nDims, int nAvgNonz, double projectionGranularity)
{
	LshTest2_Helper2<TDimIdx, TVal, TRandDirStore, TSparseVectorGenerator<TDimIdx, TVec<TKeyDat<TDimIdx, TVal> > > >(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
	LshTest2_Helper2<TDimIdx, TVal, TRandDirStore, TSparseVectorGenerator<TDimIdx, TVec<TPair<TDimIdx, TVal> > > >(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
	LshTest2_Helper2<TDimIdx, TVal, TRandDirStore, TSparseVectorGenerator<TDimIdx, THash<TDimIdx, TVal> > >(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
}

// This is like LshTest2_Helper3 + dense vectors; can only be used if TDimIdx is an integral type.
template <typename TDimIdx, typename TVal, typename TRandDirStore>
void LshTest2_Helper3a(int nKeys, int nBands, int nFuncPerBand, int nDims, int nAvgNonz, double projectionGranularity)
{
	LshTest2_Helper3<TDimIdx, TVal, TRandDirStore>(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
	LshTest2_Helper2<TDimIdx, TVal, TRandDirStore, TSparseVectorGenerator<TDimIdx, TVec<TVal> > >(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
}

// Performs a number of experiments with vectors of (TInt, TFlt) pairs and the random projection hasher;
// using various vector representations, distance measures, finalizers, and random direction stores.
int LshTest2()
{
	const int nKeys = 100, nBands = 3, nFuncPerBand = 3;
	if (false)
	{
		typedef TIntSetGenerator::TItem TItem;
		LshTest2_Helper(nKeys, TIntSetGenerator(20, 5), TMinHasher<TItem>(nBands, nFuncPerBand), TJaccardDistance<TItem>());
	}
	if (true)
	{
		int nDims = 20, nAvgNonz = 5;
		TestDots<TInt, TFlt>(nDims, nAvgNonz, 100, 3000);
		double projectionGranularity = 1;
		LshTest2_Helper3a<TInt, TFlt, TRandDirStore_Vec<TFlt> >(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
		LshTest2_Helper3a<TInt, TFlt, TRandDirStore_Hash<TInt, TFlt> >(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
		LshTest2_Helper3<TStr, TFlt, TRandDirStore_Hash<TStr, TFlt> >(nKeys, nBands, nFuncPerBand, nDims, nAvgNonz, projectionGranularity);
	}
	return 0;
}

// This class generates unit-length sparse vectors that try to simulate the statistical properties 
// of TF-IDF vectors from a real collection of text documents.
template<typename TDimIdx, typename TDat_>
class TSparseDocGenerator
{
public:
	typedef TVec<TDimIdx> TDimIdxV;
	typedef TFlt TVal;
	typedef TKeyDat<TDimIdx, TVal> TKd;
	typedef TVec<TKd> TKdV;
	typedef TDat_ TDat;
	TVec<TDimIdx> dimIndices;
	int nDims, docLen; 
	// 'alpha' is the exponent of the power-law distribution used to pick words; higher alpha means that the infrequent words are even more infrequent.
	double alpha;

	TSparseDocGenerator(int nDims_, int docLen_, double alpha_) : nDims(nDims_), docLen(docLen_), alpha(alpha_) {
		for (int dimIdx = 0; dimIdx < nDims; dimIdx++) { dimIndices.Add(); GenDimIdx(dimIdx, dimIndices.Last()); } }

protected:

	void Gen_(TRnd& rnd, TKdV& dest) const
	{
		TIntV words; words.Reserve(docLen); words.PutAll(-1);
		const double m = 1;
		const double yMax = 1 - pow(nDims / m + 1, -alpha);
		// First, we'll generate a document by selecting 'docLen' words.
		for (int i = 0; i < docLen; i++)
		{
			// We'll use a Pareto distribution, which is a continuous approximation of a Zipfian (discrete power-law) distribution.
			// http://en.wikipedia.org/wiki/Pareto_distribution
			// This distribution lives in the range [m, \inf] with the probability density p(x) = a m^a / x^{a + 1}
			// and the cumulative density P(X < x) = 1 - (m/x)^a.  Thus, if we solve y = 1 - (m/x)^a, we get
			// (1 - y)^{1/a} = m/x and x = m (1 - y)^{-1/a}.  So we can take a y that is uniform on [0, 1] and
			// obtain an x from it using this formula.  We'll round it down and subtract m to get the 0-based index of a word.
			// We want word indices in the range 0..d-1.  So if we want x - m < d, we have
			// m (1 - y)^{-1/a} - m < d, (1 - y)^{-1/a} < d/m + 1, (1 - y)^{1/a} > (d/m + 1)^{-1},
			// 1 - y > (d/m + 1)^{-a}, 1 - (d/m + 1)^{-a} > y.  This upper bound on y is curerntly in 'yMax',
			// so we'll take a y uniformly from [0, yMax] rather than from [0, 1].
			double y = rnd.GetUniDev() * yMax;
			double x = floor(m * pow(1 - y, -1.0 / alpha)) - m;
			int wordId = (int) x; if (wordId < 0) wordId = 0; if (wordId >= nDims) wordId = nDims - 1;
			words.Add(wordId);
		}
		// Now we'll convert our list of words into a TF-IDF vector.
		words.Sort();
		double norm = 0; dest.Clr();
		for (int i = 0; i < words.Len(); ) {
			// We should simulate IDFs as well.  If the probability of word x being chosen is p, then the probability
			// that it appears at least once in a document of L words is 1 - (1 - p)^L.  Thus, its document frequency
			// in a collection of N documents is going to be N * (1 - (1 - p)^L), and its IDF will be
			// log(N / DF) = log(1 / (1 - (1 - p)^L) = - log(1 - (1 - p)^L).  Now, we saw that x gets chosen if
			// x <= m (1 - y)^{- 1/a} < x + 1, so the smallest value of y at which it gets chosen is 1 - (m/x)^a and
			// the greatest is just below 1 - (m/(x + 1)^a).
			double x = words[i] + m; 
			int j = i + 1; while (j < words.Len() && words[j] == words[i]) j++;
			double tf = j - i; 
			double prob_x = pow(m / x, alpha) - pow(m / (x + 1), alpha);
			prob_x /= yMax;
			double probNotInDoc = pow(1 - prob_x, docLen);
			double df = 1 - probNotInDoc;
			double idf = (probNotInDoc < 1e-8) ? probNotInDoc : -log(df);
			double xi = tf * idf;
			norm += xi * xi;
			dest.Add(TKd(words[i], xi));
			i = j; }
		// Normalize it to unit length.
		norm = 1 / sqrt(norm); for (int i = 0; i < dest.Len(); i++) dest[i].Dat *= norm;
	}

public:

	void Init(TRnd& rnd) { }

	void Gen(TRnd& rnd, TDat& dat) const {
		TVec<TKeyDat<TInt, TVal> > kd; Gen_(rnd, kd); Convert(kd, dimIndices, dat); }

};

// This function applies the LshCollection to a set of sparse vectors generated by TSparseDocGenerator.
// It performs a number of approximate nearest-neighbor queries (using the cosine distance)
// and compares their results to the correct (i.e. non-approximate) results of the same queries.
// It also reports timings and other statistics.
int LshTest3(int nBands, int nFuncPerBand, int nResults)
{
	const int nKeys = 10000, nQueries = 1000; /*, nBands = 5, nFuncPerBand = 8, nResults = 5;*/
	//const int nKeys = 100000, nQueries = 10000; /*, nBands = 5, nFuncPerBand = 8, nResults = 5;*/
	//
	typedef TInt TDimIdx;
	typedef TFlt TVal;
	typedef TVec<TKeyDat<TDimIdx, TVal> > TDat;
	typedef TInt TKey;
	typedef TVec<TKey> TKeyV;
	typedef TCosineDistance<TDimIdx, TVal> TDistFunc;
	typedef TRandProjFinalizer_Sign<TVal, TInt> TFinalizer;
	typedef TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Set<TB32Set>, TRandDirStore_Vec<TVal>, TFinalizer> THasher;
	typedef TLshCollection<TKey, TDat, THasher, TDistFunc> TLsh;
	typedef TKeyDat<TKey, TDat> TKd;
	typedef TVec<TKd> TKdV;
	typedef TSparseDocGenerator<TDimIdx, TDat> TGen;
	typedef TDistFunc::TDistance TDistance;
	typedef TKeyDat<TDistance, TInt> TDistIntKd;
	typedef TVec<TDistIntKd> TDistIntKdV;
	//
	const int nDims = 2000, nDocLen = 200; const double alpha = 1.3;
	printf("\n\nLshTest3: %d keys, %d bands, %d functions/band; vocabulary of %d words, docLen = %d words, alpha = %g\n", nKeys, nBands, nFuncPerBand, nDims, nDocLen, alpha);
	// Generate random documents.
	TGen gen(nDims, nDocLen, alpha);
	TKdV keyDats; TRnd rnd(789);
	TTimer tm; tm.Start(); TLsh::TStDev stdNonz;
	for (int i = 0; i < nKeys; i++) {
		keyDats.Add();
		keyDats.Last().Key = i;
		gen.Gen(rnd, keyDats.Last().Dat); 
		stdNonz.Add(keyDats.Last().Dat.Len()); }
	tm.Stop(); printf("Generated %d docs in %.03f sec.  Avg nonzero components = %.2f +/- %.2f.\n", nKeys, tm.Sec(), stdNonz.Avg(), stdNonz.Std());
	// Add them to a LSH collection.
	TDistFunc distFunc;
	THasher hasher(nBands, nFuncPerBand, 0, TFinalizer());
	TLsh lsh(hasher, distFunc);
	tm.Clr(); tm.Start();
	for (int i = 0; i < nKeys; i++)
		lsh.Add(keyDats[i].Key, keyDats[i].Dat);
	tm.Stop(); printf("Added %d docs to TLsh in %.03f sec. (%.1f docs per sec)\n", nKeys, tm.Sec(), nKeys / tm.Sec());
	//
	lsh.Report(stdout);
	// Perform queries.
	TTimer tmQueries; TLsh::TStDev sumCands1, sumCands2, unionCands1, unionCands2, nBandCodes2, nAllBandCodes2, phase2Needed, totalUnionCands, stdInter, stdJaccard, stdSumRanks, stdSumRatio, stdRmse, stdDistToResults, stdIdealDist;
	TTimer tmPhase2, tmWoPhase2, tmPhase2Sort, tmPhase2Rest;
	for (int queryNo = 0; queryNo < nQueries; queryNo++)
	{
		TKeyV results; TLsh::TNnQueryReport report;
		//int iQueryDoc = rnd.GetUniDevInt(nKeys);  qDat = keyDats[iQueryDoc].Dat;
		TDat qDat; gen.Gen(rnd, qDat);
		tm.Clr(); tm.Start();
		lsh.ApproxNnQuery(qDat, nResults, 0, &results, &report);
		tm.Stop(); tmQueries.Add(tm); if (report.phase2Needed) tmPhase2.Add(tm); else tmWoPhase2.Add(tm);
		//
		sumCands1.Add(report.sumCands1);
		sumCands2.Add(report.sumCands2);
		unionCands1.Add(report.unionCands1);
		unionCands2.Add(report.unionCands2);
		nBandCodes2.Add(report.nBandCodes2);
		nAllBandCodes2.Add(report.nAllBandCodes2);
		phase2Needed.Add(report.phase2Needed ? 1 : 0);
		totalUnionCands.Add(report.unionCands1 + report.unionCands2);
		// Compute the correct results of this query.
		TDistIntKdV correctResults;
		for (int i = 0; i < nKeys; i++)
			correctResults.Add(TDistIntKd(lsh.DistFunc(qDat, keyDats[i].Dat), i));
		correctResults.Sort(); 
		TFltV rank; rank.Gen(nKeys); rank.PutAll(-1);
		for (int i = 0; i < correctResults.Len(); )
		{
			int j = i + 1; while (j < correctResults.Len() && correctResults[j].Key == correctResults[i].Key) j++;
			double avgRank = 1 + (i + (j - 1)) / 2.0;
			while (i < j) rank[correctResults[i++].Dat] = avgRank;
		}
		double sumRanks = 0, idealSum = 0, rmse = 0, distToResults = 0, idealDist = 0;
		int nInter = 0; // the intersection between our results and the correct results (if the former were also cut off at nResults)
		for (int i = 0; i < results.Len(); i++) {
			double d1 = lsh.DistFunc(qDat, keyDats[results[i]].Dat); distToResults += d1 * d1;
			double d2 = correctResults[i].Key; idealDist += d2 * d2;
			rmse += (d1 - d2) * (d1 - d2);
			double r = rank[results[i]]; IAssert(r > 0);
			if (r <= nResults) nInter++;
			sumRanks += r; idealSum += (i + 1); }
		double jaccard = nInter / float(nResults + nResults - nInter);
		stdInter.Add(nInter);
		stdJaccard.Add(jaccard);
		stdSumRanks.Add(sumRanks);
		stdSumRatio.Add(idealSum / sumRanks);
		stdRmse.Add(sqrt(rmse / nResults));
		stdDistToResults.Add(sqrt(distToResults / nResults));
		stdIdealDist.Add(sqrt(idealDist / nResults));
		if (report.phase2Needed) { tmPhase2Sort.Add(report.tmPhase2Sort); tmPhase2Rest.Add(report.tmPhase2Rest); }
	}
	printf("Performed %d queries (%d nearest neighbors) in %.03f sec (%.1f queries per sec)\n", nQueries, nResults, tmQueries.Sec(), nQueries / tmQueries.Sec());
	printf(" - queries with phase 1 only: %d queries in %.03f sec (%.1f squeries per sec)\n", int(nQueries - phase2Needed.Sum()), tmWoPhase2.Sec(), (nQueries - phase2Needed.Sum()) / tmWoPhase2.Sec());
	printf(" - queries with phase 2: %d queries in %.03f sec (%.1f squeries per sec)\n", int(phase2Needed.Sum()), tmPhase2.Sec(), phase2Needed.Sum() / tmPhase2.Sec());
	printf(" - time spent sorting in phase 2: %.03f sec total (%.2f %%)\n", tmPhase2Sort.Sec(), 100.0 * tmPhase2Sort.Sec() / tmPhase2.Sec());
	printf(" - time spent in the rest of phase 2: %.03f sec total (%.2f %%)\n", tmPhase2Rest.Sec(), 100.0 * tmPhase2Rest.Sec() / tmPhase2.Sec());
	printf(" - phase 1: sum = %.2f +/- %.2f, union = %.2f +/- %.2f\n",
		sumCands1.Avg(), sumCands1.Std(), unionCands1.Avg(), unionCands1.Std());
	printf(" - phase 2: needed in %.0f/%d queries (%.2f %%)\n", phase2Needed.Sum(), nQueries, phase2Needed.Avg() * 100);
	printf(" - phase 2: %.2f +/- %.2f lists out of %.2f +/- %.2f lists examined; sum = %.2f +/- %.2f, union = %.2f +/- %.2f, total union = %.2f +/- %.2f\n",
		nBandCodes2.Avg(), nBandCodes2.Std(),
		nAllBandCodes2.Avg(), nAllBandCodes2.Std(),
		sumCands2.Avg(), sumCands2.Std(),
		unionCands2.Avg(), unionCands2.Std(),
		totalUnionCands.Avg(), totalUnionCands.Std());
	printf(" - average intersection with the correct result set: %.2f +/- %.2f (Jaccard coefficient = %.4f +/- %.4f)\n",
		stdInter.Avg(), stdInter.Std(), stdJaccard.Avg(), stdJaccard.Std());
	printf(" - sumRanks: %.2f +/- %.2f (ideal = %.2f), ratio (bigger is better) = %.4f +/- %.4f\n",
		stdSumRanks.Avg(), stdSumRanks.Std(), 0.5 * nResults * (nResults + 1),
		stdSumRatio.Avg(), stdSumRatio.Std());
	printf(" - RMSE between the vector of distances to the query: %.4f +/- %.4f\n", stdRmse.Avg(), stdRmse.Std());
	printf(" - avg distance from q to results: %.4f +/- %.4f, ideal: %.4f +/- %.4f\n", 
		stdDistToResults.Avg(), stdDistToResults.Std(),
		stdIdealDist.Avg(), stdIdealDist.Std());
	return 0;
}

int LshTest3()
{
	LshTest3(5, 15, 5);
	return 0;
	LshTest3(5, 5, 5);
	LshTest3(5, 10, 5);
	LshTest3(5, 15, 5);
	LshTest3(10, 5, 5);
	LshTest3(10, 10, 5);
	LshTest3(10, 15, 5);
	LshTest3(15, 5, 5);
	LshTest3(15, 10, 5);
	LshTest3(15, 15, 5);
	return 0;
}

// IAsserts that the two vectors are equal.
template <typename T, typename U>
void CompareVectors(const T &t, const U &u)
{
	IAssert(t.Len() == u.Len());
	for (int i = 0; i < t.Len(); i++)
		IAssert(t[i] == u[i]);
}

// Tests TShingler::ShingleByChars on a number of random strings and different string
// representations (TStr, TChA, TUStr, char*, TIntV, etc.).  It also compares their
// results to that of ShingleByChars_Slow.
int ShingleTest()
{
	typedef TShingler<> TSh;
	TSh sh;
	TRnd rnd(543);
	TStr x1("lorem ipsum"); const TStr x2("lorem ipsum"); 
	printf("%d\n", TSh::GetLength(x1));
	printf("%d\n", TSh::GetLength(x2));
	printf("%d\n", TSh::GetLength("lorem ipsum"));
	{
		char c = (char) 0xff; uchar u = 0xff;
		int i = (int) TSh::ConvertChar(c); IAssert(i == 255);
		i = (int) TSh::ConvertChar(u); IAssert(i == 255);
	}
	const int len = 1000;
	if (!TUnicodeDef::IsDef()) { TUnicodeDef::Load("d:\\users\\janez\\dev\\NewsCluster\\NewsCluster\\dbs\\UnicodeDef.Bin"); }
	TSInt w = (short) (ushort) 0xffff;
	TSh::TInternalInt ww = TSh::ConvertChar(w);
	for (int iRange = 1; iRange < 3; iRange++)
	{
		int MaxChar = (iRange == 0) ? 255 : (iRange == 1) ? 65535 : 17 * 65536 - 1;
		//
		TIntV inStr_IntV; 
		TUIntV inStr_UIntV;
		TVec<TSInt> inStr_SIntV;
		wchar_t *inStr_wchar = new wchar_t[len + 1];
		short *inStr_short = new short[len + 1];
		char *inStr_char = new char[len + 1];
		uchar *inStr_uchar = new uchar[len + 1];
		TChA inStr_ChA;
		TStr inStr_Str;
		TUStr inStr_UStr;
		// Create several representations of the same random string.
		// iRange controls the range of character values we'll be using (which in turn determines
		// which types we may use to store the string).
		for (int i = 0; i < len; i++) { 
			int c = rnd.GetUniDevInt(1, MaxChar); IAssert(1 <= c && c <= MaxChar);
			inStr_IntV.Add(c); inStr_UIntV.Add(c); 
			if (iRange <= 1) {
				inStr_SIntV.Add((short) c);
				inStr_wchar[i] = (wchar_t) c;
				inStr_short[i] = (short) c; }
			if (iRange <= 0) {
				inStr_char[i] = (char) c;
				inStr_uchar[i] = (uchar) c;
				inStr_ChA += (char) c; }
		}
		inStr_wchar[len] = 0; inStr_short[len] = 0;
		inStr_char[len] = 0; inStr_uchar[len] = 0;
		if (iRange <= 0) inStr_Str = inStr_ChA;
		inStr_UStr = TUStr(inStr_IntV);
		// Compute shingles for various substrings of our string.  
		// We'll test that the resulting list of shingles is exactly the same regardless of the underlying representation.
		for (int ofs = 0; ofs < len; ofs++) for (int win = 1; win < len; win++)
		{
			if (win == 1) printf("%d %d    \r", iRange, ofs);
			if (ofs >= len / 5 && ofs <= 4 * len / 5 && win >= len / 5 && win <= 4 * len / 5)
				if (rnd.GetUniDevInt(10) != 7) continue;
			if (false)
			{
				TIntV vec1, vec2;
				sh.ShingleByChars(inStr_IntV, ofs, len - ofs, win, TSh::TVectorSink<TIntV>(vec1));
				sh.ShingleByChars_Slow(inStr_IntV, ofs, len - ofs, win, TSh::TVectorSink<TIntV>(vec2)); CompareVectors(vec1, vec2);
#define _(x) { vec1.Clr(); sh.ShingleByChars(inStr_ ## x, ofs, len - ofs, win, TSh::TVectorSink<TIntV>(vec1)); CompareVectors(vec1, vec2); }
				_(UIntV); 
				_(UStr);
				if (iRange <= 1) { 
					_(SIntV); 
					_(wchar); 
					_(short); }
				if (iRange <= 0) { 
					_(char); 
					_(uchar); 
					_(ChA); 
					_(Str); }
#undef _
			}
			//
			if (false) if (ofs + win <= len)
			{
				int sh1 = sh.ShingleString(inStr_IntV, ofs, win), sh2;
#define _(x) { sh2 = sh.ShingleString(inStr_ ## x, ofs, win); IAssert(sh2 == sh1); }
				_(UIntV); 
				_(UStr);
				if (iRange <= 1) { 
					_(SIntV); 
					_(wchar); 
					_(short); }
				if (iRange <= 0) { 
					_(char); 
					_(uchar); 
					_(ChA); 
					_(Str); }
#undef _
			}
		}
		//
		delete[] inStr_wchar; delete[] inStr_short; 
		delete[] inStr_char; delete[] inStr_uchar;
	}
	return 0;
}

// The following are helper classes for ShingleTest2.  Init() converts a string
// string from TIntV to a different representation (TString), whereas Done() cleans up
// any memory that may have been allocated in Init().

template<typename TChar>
struct TVecTraits
{
	typedef TVec<TChar> TString;
	static void Init(TString& s, const TIntV &src) {
		s.Gen(src.Len()); for (int i = 0; i < src.Len(); i++) s[i] = (TChar) (src[i]); }
	static void Done(TString& s) { }
};

template<typename TChar>
struct TPtrTraits
{
	typedef TChar *TString;
	static void Init(TString& s, const TIntV &src) {
		s = new TChar[src.Len() + 1]; s[src.Len()] = (TChar) 0;
		for (int i = 0; i < src.Len(); i++) s[i] = (TChar) (src[i]); }
	static void Done(TString& s) { delete[] s; }
};

struct TChATraits
{
	typedef TChA TString;
	static void Init(TString& s, const TIntV &src) {
		s.Clr(); for (int i = 0; i < src.Len(); i++) s += (char) src[i]; }
	static void Done(TString& s) { }
};

struct TStrTraits
{
	typedef TStr TString;
	static void Init(TString& s, const TIntV &src) {
		TChA buf; TChATraits::Init(buf, src); s = buf; }
	static void Done(TString& s) { }
};

struct TUStrTraits
{
	typedef TUStr TString;
	static void Init(TString& s, const TIntV &src) { s = TUStr(src); }
	static void Done(TString& s) { }
};

// Maintains a list of strings in two representations: as a TVec<TString> and a TString*.
template<typename TStringTraits_>
class TStrList
{
public:
	typedef TStringTraits_ TStringTraits;
	typedef typename TStringTraits::TString TString;
	TVec<TString> vec;
	TString *arr;

protected:
	void Put(int index, const TIntV& src) {
		TStringTraits::Init(vec[index], src);
		TStringTraits::Init(arr[index], src); }
public:
	TStrList(const TVec<TIntV> &src) { 
		int nWords = src.Len(); vec.Gen(nWords); arr = new TString[nWords];
		for (int i = 0; i < src.Len(); i++) Put(i, src[i]); }
	~TStrList() { 
		for (int i = 0; i < vec.Len(); i++) { 
			TStringTraits::Done(vec[i]);
			TStringTraits::Done(arr[i]); }
		delete[] arr; }
	// Runs ShingleByWords on both representations of the string list and triggers an assertion failure if the results differ.
	template<typename TSh>
	void Test(const TSh& sh, int ofs, int len, int win, const TIntV& correctResults) const {
		TIntV results;
		sh.ShingleByWords(vec, ofs, len, win, TSh::TVectorSink<TIntV>(results)); 
		CompareVectors(results, correctResults);
		results.Clr();
		sh.ShingleByWords(arr, ofs, len, win, TSh::TVectorSink<TIntV>(results)); 
		CompareVectors(results, correctResults); }
};

// This function tests ShingleByWords on a list of random strings in various representations
// and will trigger an assertion failure if the results do not match.
int ShingleTest2()
{
	typedef TShingler<> TSh;
	TSh sh;
	TRnd rnd(543);
	const int nWords = 1000, minWordLen = 3, maxWordLen = 10;
	if (!TUnicodeDef::IsDef()) { TUnicodeDef::Load("d:\\users\\janez\\dev\\NewsCluster\\NewsCluster\\dbs\\UnicodeDef.Bin"); }
	for (int iRange = 1; iRange < 3; iRange++)
	{
		int MaxChar = (iRange == 0) ? 255 : (iRange == 1) ? 65535 : 17 * 65536 - 1;
		// Create a few random strings.
		TVec<TIntV> strings, origStrings; strings.Gen(nWords);
		for (int i = 0; i < nWords; i++)
		{
			int len = rnd.GetUniDevInt(minWordLen, maxWordLen); 
			strings[i].Gen(len);
			for (int j = 0; j < len; j++) {
				int c = rnd.GetUniDevInt(1, MaxChar); IAssert(1 <= c && c <= MaxChar);
				strings[i][j] = c; }
		}
		origStrings = strings;
		// Create several representations of the same list of strings.
		// iRange controls the range of character values we'll be using (which in turn determines
		// which types we may use to store the string).
		TStrList<TVecTraits<TInt> > strings_IntV(strings);
		TStrList<TVecTraits<TUInt> > strings_UIntV(strings);
		TStrList<TUStrTraits> strings_UStr(strings);
		if (iRange < 1) strings.Clr();
		TStrList<TVecTraits<TSInt> > strings_SIntV(strings);
		TStrList<TPtrTraits<wchar_t> > strings_wchar(strings);
		TStrList<TPtrTraits<short> > strings_short(strings);
		if (iRange < 2) strings.Clr();
		TStrList<TPtrTraits<char> > strings_char(strings);
		TStrList<TPtrTraits<uchar> > strings_uchar(strings);
		TStrList<TChATraits> strings_ChA(strings);
		TStrList<TStrTraits> strings_Str(strings);
		// Compute shingles for various substrings of our string.  
		// We'll test that the resulting list of shingles is exactly the same regardless of the underlying representation.
		for (int ofs = 0; ofs < nWords; ofs++) for (int win = 1; win < nWords; win++)
		{
			if (ofs >= nWords / 5 && ofs <= 4 * nWords / 5 && win >= nWords / 5 && win <= 4 * nWords / 5)
				if (rnd.GetUniDevInt(10) != 7) continue;
			TIntV vec1;
			sh.ShingleByWords(origStrings, ofs, nWords - ofs, win, TSh::TVectorSink<TIntV>(vec1));
#define _(x) (strings_ ## x).Test(sh, ofs, nWords - ofs, win, vec1)
			_(IntV); 
			_(UIntV); 
			_(UStr);
			if (iRange <= 1) { 
				_(SIntV); 
				_(wchar); 
				_(short); }
			if (iRange <= 0) { 
				_(char); 
				_(uchar); 
				_(ChA); 
				_(Str); }
		}
	}
	return 0;
}

// A minimalistic example of how to use the TLshCollection while getting as little involved
// as possible with its internals.
//
// In this example, our documents will be described by a TIntFltKdV feature vector, and each
// document is uniquely identified by a key, which will be a TStr.  The cosine distance is
// used for nearest-neighbor queries.
void LshExample()
{
	int nBands = 5, nFuncPerBand = 5;
	Lsh_IntFltKdV_Cosine_8<TStr>::PLsh lsh = Lsh_IntFltKdV_Cosine_8<TStr>::New(nBands, nFuncPerBand);
	TIntFltKdV v; // v = ... (initialization of v omitted in this example)
	lsh->Add("doc1", v);
	// v = ... (omitted in this example)
	lsh->Add("doc2", v);
	int nResults = 10; TStrV resultKeys; 
	lsh->ApproxNnQuery(v, nResults, resultKeys);
	lsh->Del("doc1");
}

// Another example of using TLshCollection.  In this case, each document has a unique key (TStr)
// and the full text (also TStr); the full text will be shingled into a TIntV, which can then
// be used in the TLshCollection.  
void LshExample2()
{
	int nBands = 5, nFuncPerBand = 5, shingleLen = 4;
	Lsh_IntSet_Jaccard<TStr>::PLsh lsh = Lsh_IntSet_Jaccard<TStr>::New(nBands, nFuncPerBand);
	TStr doc1text, doc2text, queryText; // initialization omitted in this example
	TShingler<> shingler; TIntV v; 
	shingler.ShingleByChars(doc1text, shingleLen, v); lsh->Add("doc1", v);
	shingler.ShingleByChars(doc2text, shingleLen, v); lsh->Add("doc2", v);
	int nResults = 10; TStrV resultKeys;
	shingler.ShingleByChars(queryText, shingleLen, v); 
	lsh->ApproxNnQuery(v, nResults, resultKeys);
	lsh->Del("doc1");
}



} // LshTest