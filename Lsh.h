#ifndef ____LSH_H_INCLUDED____
#define ____LSH_H_INCLUDED____

// Locality sensitive hashing
// Based on ch. 3 of Rajaraman, Leskovec, Ullman, "Mining of Massive Datasets"

//-----------------------------------------------------------------------------
// TLshCollection
//-----------------------------------------------------------------------------
//
// This class maintains a set of (key, dat) pairs, where 'key' is a unique
// identifier of each pair.  The main feature of the class is that it uses
// locality-sensitive hashing to support (approximate) nearest-neighbor queries 
// based on the dat; i.e. given a TDat q, you can find keys whose corresponding 
// dats are close to q.
// 
// The hashing itself is performed by THasher, which must be capable of 
// calculating a b-tuple of "band codes" from any given TDat (see the description
// of IHasher and the various concrete hasher templates below).
// An object x is considered as a candidate near neighbor of q if they match
// in at least one of the band codes.  To make the detection of such matches
// easier, TLshCollection maintains b hash tables (TBandHash) which map
// a band code into a chain (i.e. a doubly linked list) of key IDs of all objects 
// in the collection with this particular band code in this particular band.
//
// The keyHash hash table maps each key into a TKeyRec structure, which contains
// the dat corresponding to this key as well as a vector of TKeyInBand structures.
// These contain, for each band, the location of the key in that band's bandHash 
// table and the previous/next key ID in the doubly-linked list of keys with the same 
// band code for that particular band.
//
// For each band, we have a bandHash table which maps band codes into TBandCodeRec
// structures.  Such a structure simply gives the key ID of the first and last key
// with this particular band code in this particular band.  The rest of the keys
// with this band code can be reached by following the 'prev' and 'next' fields
// in the keyHash table.
//
// TLshCollection supports:
// - Add(key, dat) - returns keyId
// - Del(key), DelKeyId(keyId)
// - ApproxNnQuery(qDat, nResults, dest) - returns the 'nResults' keys whose
//   corresponding dat's are (approximately) closest to 'qDat', in increasing 
//   order of distance from qDat.  
// - ApproxNnQuery_KeyIds(qDat, nResults, dest) - same as above, but returns
//   key IDs (into keyHash) instead of keys.
//
// Typical scenarios:
//
// (1) Sets of integers:
//     TDat = TIntV, containing the elements of a set of integers in increasing order
//     TDistFunc = TJaccardDistance
//     THasher = TMinHasher
// 
// (2) Sparse vectors (e.g. TFIDF vectors), cosine similarity:
//     TDat = TIntFltKdV, TIntFltPrV, TIntFltH, even a dense vector TFltV
//            or TStrFltKdV, TStrFltPrV, TStrFltH
//     TDistFunc = TCosineDistance
//     THasher = TRandProjHasher with TBandCodeAccesss_Set and TRandProjFinalizer_Sign
// 
// (3) Sparse vectors, Euclidean distance:
//     TDat = same as (2)
//     TDistFunc = TEuclideanDistance
//     THasher = TRandProjHasher with TBandCodeAccess_Vec and TRandProjFinalizer_Floor
// 
// See the classes at the end of this file for typedef's that cover these typical scenarios:
// Lsh_IntSet_Jaccard, Lsh_IntFltKdV_Cosine_8 (and _32), Lsh_IntFltKdV_Euclidean.
//    
//
// A few notes about the choice of parameters:
//
// - The number of bands (b): when performing an approximate-NN query, we examine
//   one chain of candidates from each band (namely the chain consisting of those
//   object that match the query object's band code in this particular band).
//   So the more bands there are, the longer it will take to perform the
//   NN query, but the results are also likely to be better (the more candidates
//   we examine, the closer the neighbors we are likely to find).
//
// - Each band code is usually itself an r-tuple of hash codes computed by
//   r different hash functions.  The LshCollection doesn't actually know about
//   these internals of a band code, but the THasher does.  If we increase r,
//   the number of distinct possible band codes also increases (exponentially),
//   and so the average length of a chain in each bandHash table (i.e. the list
//   of objects with the same band code in this particular band) decreases
//   accordingly.  If the chains are too long, the NN queries will be slower
//   (since the query has to examine one chain from each bandHash table),
//   but if the chains are too short, the results of the NN queries will be poor
//   (because too few candidate near neighbors will be examined).  Thus, the
//   choice of r should try to ensure that we fall somewhere between these two extremes.
//
// Here is a concrete example from a synthetic dataset of 100000 objects
// (TIntFltKdV document vectors with an average of 16.6 nonzero components
// from a 2000-dimensional space):
//
// b (# of    r (# of      approx.     avg       queries                  avg dist
// bands)    hashes per    # of       chain      per sec    unionCands1   to query
//             band)       chains    length      
// -------   ----------    --------  --------    --------   -----------   --------
//    5          5             32     3125          40.5      20077        0.2181
//    5         10           1024       97.6        76.8       1146        0.2636
//    5         15          22000        4.5       284.4         73.3      0.3792
//   10          5             32     3125          20.9      36409        0.2116
//   10         10           1024       97.6       171.3       2020        0.2346
//   10         15          22000        4.5       298.4        143.4      0.3075
//   15          5             32     3125          15.2      47610        0.2106
//   15         10           1024       97.6       138.8       3157        0.2238
//   15         15          22000        4.5       543.3        203.0      0.2786
//
// Notes:
// - b (# of bands) = the number of bands, as returned by THasher::GetNumBands()
// - r (# of hashes per band) = this means that each band code was an r-tuple of hash codes.
//   In this experiment, the hash codes were of the form h(x) = sign(w^t x) for various
//   random directions w; so with r hash codes and 2 possible values for each hash code,
//   this means there are up to 2^r possible distinct band codes.
// - approx # of chains: this is the number of distinct band codes that actually occur
//   in the dataset.  This number can of course vary from one band to another;
//   in the case of r = 5 and 10, we actually always saw all possible 2^r band codes.
// - avg chain length: this is basically the number of objects (100000) divided by
//   the number of chains (from the previous column);
// - unionCands1: the approximate NN search examines one chain from each bandHash
//   table, and unionCands1 is the total number of distinct candidates examined in this way.
//   This is basically (# of bands) * (avg chain length), or actually a bit less than this
//   because the same candidate might appear in several of these chains (if it matches
//   our query point in several band codes).    [For more details, see the description
//   of unionCands1 in TNnQueryReport below.]
// - queries per sec: how many nearest-neighbor queries (asking for 5 nearest neighbors)
//   were performed, on average, per second.  The time spent for each query is roughly
//   linearly proportional to the number of candidates examined.
// - avg dist to query: for each query, we compute the average cosine distance to the
//   5 neighbours returned by our approximate NN query method, and finally compute the
//   average of this over all 10000 queries.  If we performed an accurate (i.e.
//   non-approximate) nearest neighbor query, the resulting average would be 0.2101.
//   The closer we are to this, the better our approximate NN query results are.
//
// As we can see from this, we have to choose b and r in such a way as to find a suitable
// balance.  If we cause the NN query to examine too few candidates, results will be worse,
// but if we cause it to examine too many candidates, the queries will take too much time.
//
// Another thing that has a huge impact on the speed of NN queries is whether
// the second phase of the query process is required.  In the first phase, we examine
// those candidates that match our query object in at least one band code.  As we saw
// above, there are unionCands1 such candidates; if this is less than nResults (the number
// of neighbors requested by the user -- in our experiment here, we had nResults = 5),
// the query method then proceeds to phase 2, in which it examines more and more chains
// (which match our query object's band codes only approximately) until it accumulates
// enough neighbors.  PHASE 2 IS SIGNIFICANTLY SLOWER THAN PHASE 1, and should typically
// be avoided.  If you find that the queries often proceed to phase 2, this might be a
// sign that the chains are too short and you should decrease r.  
// - You can pass a TNnQueryReport instance to ApproxNnQuery to receive statistics
//   about the number of candidates examined and whether phase 2 was needed or not.
// - You can also tell ApproxNnQuery to ignore phase 2 altogether and return the results 
//   from phase 1 even if there are fewer than nResults of them.
//
// In the above experiment, phase 2 was needed only in 28 out of 10000 queries at b = 5, r = 15
// and in 1 out of 10000 queries at b = 10, r = 15.  

template<
	typename TKey_, // should be a unique identifier of each object
	typename TDat_, // this is the object used in measuring distances etc. - could be a vector, a set of strings, etc.
	typename THasher_, // to see what we expect from the hasher, see the comments on the IHasher template below 
	typename TDistFunc_, // measure of distance between TDat's
	typename TKeyHashFunc_ = TDefaultHashFunc<TKey_> // used to put keys into a hash table (keyHash)
>
class TLshCollection
{
public:
	TCRef CRef;
private:
	
	// Instantiating this template causes a compile-time error if T and U aren't the same class.
	template<typename T, typename U> class TSameClass { private: TSameClass(); };
	template<typename T> class TSameClass<T, T> { public: TSameClass() { } };

public:

	typedef TKey_ TKey;
	typedef TDat_ TDat;
	typedef TVec<TKey> TKeyV;
	typedef THasher_ THasher;
	typedef TDistFunc_ TDistFunc;
	typedef typename THasher::TBandCode TBandCode; // usually a vector of hash codes from all the functions in the band; but if each hash code returns 1 bit, this might be just an int
	typedef TVec<TBandCode> TBandCodeV;
	typedef typename TDistFunc::TDistance TDistance; // might be int or float etc.
	typedef typename THasher::TBandCodeHashFunc TBandCodeHashFunc;
	typedef TKeyHashFunc_ TKeyHashFunc;

	// For each band code of each band, we have a doubly-linked list of TKeyInBand structures
	// for all the keys that were hashed into this particular band code for this particular band.
	// The beginning and end of the list are pointed to by a TBandCodeRec structure (see below).
	struct TKeyInBand
	{
		int bandCodeId; // keyId in the band hash where this object has been hashed into
		int prev, next; // keyIds into 'keyHash' for the previous/next object with the same band code

		TKeyInBand() : bandCodeId(-1), prev(-1), next(-1) { }
		TKeyInBand(TSIn& SIn) { SIn.Load(bandCodeId); SIn.Load(prev); SIn.Load(next); }
		void Save(TSOut& SOut) const { SOut.Save(bandCodeId); SOut.Save(prev); SOut.Save(next); }
	};

	typedef TVec<TKeyInBand> TKeyInBandV;

	struct TKeyRec
	{
		TDat dat;
		TKeyInBandV bands; // indexed by band no.

		TKeyRec() { }
		TKeyRec(TSIn& SIn) : dat(SIn), bands(SIn) { }
		void Save(TSOut& SOut) const { dat.Save(SOut); bands.Save(SOut); }
	};

	struct TBandCodeRec
	{
		int first, last; // keyIds into 'keyHash' for the first/last object with the same band code
		TBandCodeRec() : first(-1), last(-1) { }
		TBandCodeRec(TSIn& SIn) { SIn.Load(first); SIn.Load(last); }
		void Save(TSOut& SOut) const { SOut.Save(first); SOut.Save(last); }
	};

	typedef TPair<TDat, TIntV> TRecord;
	typedef THash<TKey, TKeyRec, TKeyHashFunc> TKeyRecH;
	typedef THash<TBandCode, TBandCodeRec, TBandCodeHashFunc> TBandHash; // key = a band code; keyId = will be referred to as a 'band code ID' to reduce confusion; dat = a linked list of objects whose bandCode for this particular band was 'key'
	typedef TVec<TBandHash> TBandHashV;

	
protected:

	THasher hasher;
	TDistFunc distFunc;
	TKeyRecH keyHash; // key: TKey; dat: pair (TDat, TIntV), where the IntV contains keyId's of the records referring to this key in the individual bandHashes.  This makes deletion easier.
	TBandHashV bandHashes; 

public:

	void Clr() {
		keyHash.Clr();
		bandHashes.Clr();
		bandHashes.Gen(hasher.GetNumBands()); }

	void Clr(const THasher& newHasher) {
		hasher = newHasher; Clr(); }

	TLshCollection(const THasher& Hasher, const TDistFunc &DistFunc) : hasher(Hasher), distFunc(DistFunc) { Clr(); }
	TLshCollection(const THasher& Hasher) : hasher(Hasher) { Clr(); }
	TLshCollection(TSIn& SIn) : hasher(SIn), keyHash(SIn), bandHashes(SIn) { SIn.LoadCs(); }
	void Save(TSOut& SOut) const { hasher.Save(SOut); keyHash.Save(SOut); bandHashes.Save(SOut); SOut.SaveCs(); }

	bool Empty() const { return keyHash.Empty(); }
	int Len() const { return keyHash.Len(); }
	int GetKeyId(const TKey& key) const { return keyHash.GetKeyId(key); }
	bool IsKey(const TKey& key) const { return keyHash.IsKey(key); }
	const TKey& GetKey(int keyId) const { return keyHash.GetKey(keyId); }
	// Note: you should *NOT* modify the TDat returned by GetDat, at least not in such
	// a way as would cause its band codes to change.  If you need to modify a dat, simply 
	// call Add again -- this will ensure that the band hash tables get updated properly.
	TDat& GetDat(const TKey& key) { return keyHash.GetDat(key).dat; }
	const TDat& GetDat(const TKey& key) const { return keyHash.GetDat(key).dat; }

	// These methods provide read-only access to the keyHash table, the hasher, and the distance function.
	// The user of the TLshCollection isn't really meant to access these things, and certainly not to modify them.
	const TKeyRecH &KeyHash() const { return keyHash; }
	const THasher &Hasher() const { return hasher; }
	const TDistFunc &DistFunc() const { return distFunc; }
	TDistance DistFunc(const TDat& dat1, const TDat& dat2) const { return distFunc(dat1, dat2); }

protected:

	// Deletes entries referring to the given keyId from the band hashes.
	void DelFromBandHashes(const int keyId)
	{
		const int nBands = hasher.GetNumBands(); 
		const TKeyInBandV &kibv = keyHash[keyId].bands;
		Assert(kibv.Len() == nBands);
		for (int iBand = 0; iBand < nBands; iBand++) {
			const TKeyInBand &kib = kibv[iBand];
			const int bandCodeId = kib.bandCodeId, prev = kib.prev, next = kib.next;
			// For the iBand'th band, our object (which appears in keyHash under the ID 'keyId') has the band code
			// bandHashes[iBand].GetKey(bandCodeId); the corresponding dat there is a linked list of key IDs and
			// we have to remove our current keyId from that linked list.
			TBandCodeRec &bcRec = bandHashes[iBand][bandCodeId];
			if (prev < 0) { 
				IAssert(bcRec.first == keyId); 
				bcRec.first = next; }
			else { 
				Assert(keyHash[prev].bands[iBand].next == keyId); 
				Assert(keyHash[prev].bands[iBand].bandCodeId == bandCodeId); 
				keyHash[prev].bands[iBand].next = next; }
			if (next < 0) { 
				IAssert(bcRec.last == keyId); 
				bcRec.last = prev; }
			else { 
				Assert(keyHash[next].bands[iBand].prev == keyId); 
				Assert(keyHash[next].bands[iBand].bandCodeId == bandCodeId); 
				keyHash[next].bands[iBand].prev = prev; }}
	}

public:

	// Adds a key with the corresponding dat.  If the key was already in the keyHash, it will remain
	// where it was, but its corresponding dat will be modified (and the band hashes will be suitably updated).
	int Add(const TKey& key, const TDat& dat)
	{
		const int nBands = hasher.GetNumBands();
		// If this key is already in the keyHash, it is also in the band hashes, in which case
		// it needs to be removed from there (since we have a new dat now and its position in the band
		// hashes will therefore likely change).  If the key isn't in the keyHash yet, we'll add it now.
		int keyId = keyHash.GetKeyId(key);
		if (keyId >= 0) 
			DelFromBandHashes(keyId);
		else {
			keyId = keyHash.AddKey(key);
			keyHash[keyId].bands.Gen(nBands); }
		keyHash[keyId].dat = dat;
		TKeyInBandV &kibv = keyHash[keyId].bands;
		Assert(kibv.Len() == nBands);
		kibv.PutAll(TKeyInBand()); 
		// Now compute its band codes and place it into the suitable linked lists in the band hashes.
		TBandCodeV bandCodes; 
		hasher.CalcBandCodes(dat, bandCodes);
		for (int iBand = 0; iBand < nBands; iBand++)
		{
			TBandHash &bandHash = bandHashes[iBand];
			int bandCodeId = bandHash.AddKey(bandCodes[iBand]);
			// For the iBand'th band, our object has the band code 'bandCodes[iBand]', whose ID in bandHashes[iBand] is 'bandCodeId'.
			// Add keyId into the corresponding linked list.
			TBandCodeRec &bcRec = bandHash[bandCodeId];
			TKeyInBand &kib = kibv[iBand];
			kib.bandCodeId = bandCodeId; kib.prev = -1; kib.next = bcRec.first;
			if (bcRec.first >= 0) { 
				Assert(keyHash[bcRec.first].bands[iBand].bandCodeId == bandCodeId); 
				Assert(keyHash[bcRec.first].bands[iBand].prev < 0); 
				keyHash[bcRec.first].bands[iBand].prev = keyId; }
			bcRec.first = keyId; if (bcRec.last < 0) bcRec.last = keyId;
		}
		return keyId;
	}

	void DelKeyId(const int keyId) {
		Assert(keyHash.IsKeyId(keyId));
		DelFromBandHashes(keyId);
		keyHash.DelKeyId(keyId); }

	// If the key exists, Del deletes it and return true; otherwise it returns false.
	bool Del(const TKey& key) {
		int keyId = keyHash.GetKeyId(key);
		if (keyId < 0) return false;
		DelKeyId(keyId); return true; }

protected:

	// The heap ensures that cmp(p, c) is true for each (parent, child) pair in the heap.
	template<typename TKey_, typename TDat_, typename TCmp_ = TLss<TKey_> >
	class THeap
	{
	public:
		typedef TKey_ TKey; typedef TDat_ TDat; typedef TCmp_ TCmp;
		// Note: we use TPairs instead of TKeyDat's for the following reason.  At the end of
		// an approximate nearest-neighbor query, we will call v.Sort() to sort the results in
		// ascending order of distance from the query point.  If several results are at the same distance,
		// then having TPairs will cause ties to be broken by sorting on the second component of the pair,
		// which is the keyId of each result.  Having TKeyDat's would mean that the entries with
		// the same distance are returned in a random order, and what is worse, this order can change
		// from one call to the next (because we sort by calling v.Sort() and this ultimately uses
		// a random number generator from a global variable, TInt::Rnd).
		typedef TPair<TKey, TDat> TEntry;
		TVec<TEntry> v;
		TCmp cmp;
		void Clr() { v.Clr(); }
		int Len() const { return v.Len(); }
		int Sift(int i) {
			int n = v.Len(); TEntry ent = v[i];
			while (2 * i + 1 < n) {
				int ci = 2 * i + 1;
				if (ci + 1 < n && cmp(v[ci + 1].Val1, v[ci].Val1)) ci += 1;
				if (! cmp(v[ci].Val1, ent.Val1)) break;
				v[i] = v[ci]; i = ci; }
			v[i] = ent; return i; }
		int Lift(int i) {
			TEntry ent = v[i];
			while (i > 0) {
				int p = (i - 1) / 2;
				if (! cmp(ent.Val1, v[p].Val1)) break;
				v[i] = v[p]; i = p; }
			v[i] = ent; return i; }
		int Add(const TKey& key, const TDat& dat) { int i = v.Add(); v[i].Val1 = key; v[i].Val2 = dat; return Lift(i); }
		bool Empty() const { return v.Empty(); }
		void DelRoot() { IAssert(! v.Empty()); v[0] = v.Last(); v.DelLast(); if (! v.Empty()) Sift(0); }
	};

	typedef THeap<TDistance, TInt /* keyId */, TGtr<TDistance> > TNnHeap;

	// Computes the distance from 'keyId' to the query object 'dat' and adds it to the heap
	// (except if 'maxResults' better results are already in the heap).
	void NnQuery_ProcessCandidate(const TDat& dat, int maxResults, const int keyId, TNnHeap &heap) const
	{
		TDistance dist = distFunc(dat, keyHash[keyId].dat);
		if (heap.Len() < maxResults || heap.cmp(heap.v[0].Val1, dist)) {
			if (heap.Len() >= maxResults) heap.DelRoot();
			heap.Add(dist, keyId); }
	}

public:

	// An utility class to compute averages and standard deviations more easily.
	struct TStDev
	{
		int n; double sum, sum2, min_, max_;
		TStDev() : n(0), sum(0), sum2(0), min_(0), max_(0) { }
		void Add(double x) { 
			if (n == 0) min_ = x, max_ = x;
			n += 1; sum += x; sum2 += x * x;
			if (x < min_) min_ = x; 
			if (x > max_) max_ = x; }
		int Count() const { return n; }
		double Sum() const { return sum; }
		double Avg() const { if (n <= 0) return 0.0; else return sum / double(n); }
		double Var() const { if (n <= 1) return 0.0;
			double avg = sum / double(n); double s2 = sum2 / double(n) - avg * avg;
			return (s2 < 1e-8) ? 0.0 : s2; }
		double Std() const { if (n <= 1) return 0.0;
			double avg = sum / double(n); double s2 = sum2 / double(n) - avg * avg;
			return (s2 < 1e-8) ? 0.0 : sqrt(s2); }
		double Min() const { return min_; }
		double Max() const { return max_; }
	};

	// One of the versions of ApproxNnQuery can return this structure with potentially
	// interesting diagnostic information about the amount of effort involved in processing the query.
	struct TNnQueryReport
	{
		// The approximate NN query starts by computing the band codes of the query object
		// and examining the lists of objects that match the query object in at least one band code.
		// 'sumCands1' receives the sum of the lengths of those lists, 'unionCands1' receives
		// the number of objects in the union of those lists (this is also the number of calls
		// made to NnQuery_ProcessCandidate).
		int sumCands1, unionCands1;
		// If 'unionCands1' was less than the number of nearest neighbors requested by
		// the user, we'll have to continue with phase 2 (looking through lists of candidates
		// that match the query object's band codes only approximately, not entirely).
		bool phase2Needed;
		// In phase 2, we have to examine all band codes in all the band hash tables 
		// and sort them in decreasing order of similarity to the query object's band codes.
		int nAllBandCodes2;
		// Then we keep examining these band codes in this order, and for each band code
		// we examine the corresponding list of objects.  'nBandCodes2' receives the number
		// of these band codes, 'sumCands2' receives the sum of the lengths of those lists,
		// and 'unionCands2' receives the size of the union of those lists, excluding any
		// candidates that were already processed in phase 1.
		int nBandCodes2, sumCands2, unionCands2;
#ifdef TTimer_DEFINED
		TTimer tmPhase2Sort, tmPhase2Rest;
#endif

		void Clr() { sumCands1 = 0; unionCands1 = 0; phase2Needed = false; nAllBandCodes2 = 0; nBandCodes2 = 0; sumCands2 = 0; unionCands2 = 0; 
#ifdef TTimer_DEFINED
			tmPhase2Sort.Clr(); tmPhase2Rest.Clr(); 
#endif
		}
		TNnQueryReport() { Clr(); }
		void Print(FILE *f = 0) const {
			if (! f) f = stdout;
			fprintf(f, "phase 1: sum = %d, union = %d; phase2: ", sumCands1, unionCands1);
			if (! phase2Needed) fprintf(f, "not needed");
			else fprintf(f, "%d/%d lists examined, sum = %d, union = %d, total union = %d",
				nBandCodes2, nAllBandCodes2, sumCands2, unionCands2, unionCands1 + unionCands2); }
	};

public:

	// This method performs an approximate nearest neighbor query, looking for 'nResults' nearest
	// neighbors of the object 'dat'.  They are returned, in increasing order of distance from 'dat',
	// in '*destKeyIds' and '*destKeys'.  Old contents of these vectors are cleared before the new
	// results are added.  One or both of these pointers may be null if the caller doesn't
	// need key IDs and/or keys.  
	//
	// '*report' receives some information useful for profiling (see the TNnQueryReport class 
	// for more details).  It can be set to 0 if the user doesn't need this information.
	//
	// Normally the function uses all the band hash tables to find candidate near neighbors, 
	// but you can use 'nBands_' to limit it to the first few band hash tables to save time 
	// (at the risk of missing some potentially good candidates).  
	//
	// The nearest-neighbor search consists of two phases.  In phase 1 we look for candidates 
	// that match 'qDat' in at least one band code.  If less than 'nResults' candidates are found,
	// we also perform phase 2, in which we look for candidates that match one of 'qDat's band
	// codes at least approximately.  This phase can be a lot more time-consuming, so you can 
	// set 'skipPhase2' to true to skip it (in this case the method might return less than 
	// 'nResults' results even if there are more than 'nResults' objects in the collection).
	//
	// The return value is the number of results returned.  
	//
	// Note: one might imagine that a function like ApproxNnQuery should be marked as const;
	// but it needs to call hasher.CalcBandCodes, which cannot be const since some hashers
	// can modify themselves while calculating hash codes (e.g. the random projection hasher,
	// if it needs to add a new dimension that it hasn't encountered before).
	//
	int ApproxNnQuery(const TDat& dat, int nResults, TIntV *destKeyIds, TKeyV *destKeys, TNnQueryReport *report, int nBands_ = -1, bool skipPhase2 = false) 
	{
		int nBands = hasher.GetNumBands();
		if (nBands_ > 0 && nBands_ < nBands) nBands = nBands_;
		TBandCodeV bandCodes; hasher.CalcBandCodes(dat, bandCodes);
		if (report) report->Clr();
		// First, let's examine objects that match with 'dat' completely in at least one band code.
		THash<TInt, TVoid> cands; 
		TNnHeap heap;
		for (int iBand = 0; iBand < nBands; iBand++)
		{
			const TBandHash &bandHash = bandHashes[iBand];
			int bandCodeId = bandHash.GetKeyId(bandCodes[iBand]);
			if (bandCodeId < 0) continue;
			for (int keyId = bandHash[bandCodeId].first; keyId >= 0; ) {
				cands.AddKey(keyId);
				if (report) report->sumCands1++;
				Assert(keyHash[keyId].bands[iBand].bandCodeId == bandCodeId);
				keyId = keyHash[keyId].bands[iBand].next; }
		}
		if (report) report->unionCands1 = cands.Len();
		for (int i = cands.FFirstKeyId(); cands.FNextKeyId(i); ) 
			NnQuery_ProcessCandidate(dat, nResults, cands.GetKey(i), heap);
		// The idea of locality-sensitive hashing is that this should have provided us with enough candidates to
		// find a good approximate answer to our nearest-neighbor query.  If we didn't find enough candidates,
		// we can start looking through objects that don't match with 'dat' completely in any band code, but
		// that at least match it approximately.  However, this can quickly become much more inefficient.
		if (heap.v.Len() < nResults && ! skipPhase2)
		{
			if (report) report->phase2Needed = true;
			typedef TIntPr TBandNoBandCodePr;
			// The following is a pair and not a keydat to ensure that the results are completely repeatable.
			// If you use a keydat, the order in which entries with the same key are arranged by TVec::Sort
			// is basically random and unpredictable (because Partition uses a global random generator, TInt::Rnd).
			typedef TPair<TInt, TBandNoBandCodePr> TDistAndBandCodeKd;
			TVec<TDistAndBandCodeKd> v;
#ifdef TTimer_DEFINED
			if (report) report->tmPhase2Sort.Start();
#endif
			for (int iBand = 0; iBand < nBands; iBand++)
			{
				const TBandHash &bandHash = bandHashes[iBand];
				if (report) report->nAllBandCodes2 += bandHash.Len();
				for (int bandCodeId = bandHash.FFirstKeyId(); bandHash.FNextKeyId(bandCodeId); )
				{
					int dist = hasher.BandCodeDistance(bandCodes[iBand], bandHash.GetKey(bandCodeId));
					Assert(dist >= 0); if (dist == 0) continue;
					v.Add(TDistAndBandCodeKd(dist, TBandNoBandCodePr(iBand, bandCodeId)));
				}
			}
			v.Sort();
#ifdef TTimer_DEFINED
			if (report) report->tmPhase2Sort.Stop();
			if (report) report->tmPhase2Rest.Start();
#endif
			for (int i = 0; i < v.Len(); i++)
			{
				int iBand = v[i].Val2.Val1, bandCodeId = v[i].Val2.Val2;
				if (report) report->nBandCodes2++;
				for (int keyId = bandHashes[iBand][bandCodeId].first; keyId >= 0; keyId = keyHash[keyId].bands[iBand].next) {
					Assert(keyHash[keyId].bands[iBand].bandCodeId == bandCodeId);
					if (report) report->sumCands2++;
					if (cands.IsKey(keyId)) continue;
					cands.AddKey(keyId); // Note: this prevents us from examining the same candidate multiple times (and possibly
						// risk having the same key included multiple times in 'heap' if the distance to 'dat' is not exactly
						// the same each time (though it should be), but the downside is that the 'cands' hash table could potentially
						// grow to include all keys in the collection, so there would be a high memory consumption.
						// An alternative would be to maintain a hash table of all the keys in 'heap' to make sure that no
						// key is included in 'heap' multiple times, but we would still potentially have to consider the same
						// key multiple times (if we find it in different bandHashes).
					NnQuery_ProcessCandidate(dat, nResults, keyId, heap); 
					if (report) report->unionCands2++; }
				if (heap.Len() >= nResults && (i == v.Len() - 1 || v[i].Val1 < v[i + 1].Val1))
					break; // we've seen enough 
			}
#ifdef TTimer_DEFINED
			if (report) report->tmPhase2Rest.Stop();
#endif
		}
		if (false) { printf("\n Heap:"); for (int i = 0; i < heap.v.Len(); i++) printf(" %d", int(heap.v[i].Val2)); printf(" -> "); }
		// Return the results;
		int n = heap.v.Len();
		heap.v.Sort();
		if (false) { for (int i = 0; i < heap.v.Len(); i++) printf(" %d", int(heap.v[i].Val2)); printf("\n"); }
		if (destKeyIds) { destKeyIds->Clr(); destKeyIds->Gen(n); }
		if (destKeys) { destKeys->Clr(); destKeys->Gen(n); }
		for (int i = 0; i < heap.v.Len(); i++) {
			int keyId = heap.v[i].Val2;
			if (destKeyIds) (*destKeyIds)[i] = keyId;
			if (destKeys) (*destKeys)[i] = keyHash.GetKey(keyId); }
		return n;
	}

	// Wrappers around the generalized version of ApproxNnQuery.
	int ApproxNnQuery(const TDat& dat, int nResults, TKeyV &dest) { 
		return ApproxNnQuery(dat, nResults, 0, &dest, 0); }
	int ApproxNnQuery_KeyIds(const TDat& dat, int nResults, TIntV &dest) {
		return ApproxNnQuery(dat, nResults, &dest, 0, 0); }

	// Makes a bunch of paranoid IAssert'ions regarding the contents of our data structures.
	void Validate() 
	{
		const int nBands = hasher.GetNumBands();
		const int nKeys = keyHash.Len();
		for (int keyId = keyHash.FFirstKeyId(); keyHash.FNextKeyId(keyId); )
		{
			const TKeyRec &rec = keyHash[keyId];
			IAssert(rec.bands.Len() == nBands);
			// Recalculate the band codes for this key, look them up in the bandHashes and
			// and make sure they match the ones that are currently stored in this key's entry
			// in the keyHash table.
			TBandCodeV bandCodes; hasher.CalcBandCodes(rec.dat, bandCodes);
			IAssert(bandCodes.Len() == nBands);
			for (int iBand = 0; iBand < nBands; iBand++)
			{
				int bandCodeId = bandHashes[iBand].GetKeyId(bandCodes[iBand]);
				IAssert(bandCodeId >= 0);
				IAssert(rec.bands[iBand].bandCodeId == bandCodeId);
				const TBandCode &bandCode = bandHashes[iBand].GetKey(bandCodeId);
				IAssert(bandCode == bandCodes[iBand]);
				IAssert(hasher.BandCodeDistance(bandCode, bandCodes[iBand]) == 0);
			}
		}
		// Each bandHash consists of a bunch of (bandCode, linked list) pairs. 
		// Traverse all the lists and make sure they are structurally sound, and that
		// each key appears in exactly one list for each band.
		for (int iBand = 0; iBand < nBands; iBand++)
		{
			THash<TInt, TVoid> keySeen;
			const TBandHash &bandHash = bandHashes[iBand];
			for (int bandCodeId = bandHash.FFirstKeyId(); bandHash.FNextKeyId(bandCodeId); )
			{
				TBandCode bandCode = bandHash.GetKey(bandCodeId);
				const TBandCodeRec &bcr = bandHash[bandCodeId];
				int prev = -1, keyId = bcr.first;
				while (keyId >= 0) {
					IAssert(! keySeen.IsKey(keyId));
					keySeen.AddKey(keyId);
					IAssert(keyHash.IsKeyId(keyId));
					const TKeyInBand &kib = keyHash[keyId].bands[iBand];
					IAssert(kib.bandCodeId == bandCodeId);
					IAssert(kib.prev == prev);
					prev = keyId; keyId = kib.next; }
				IAssert(prev == bcr.last);
			}
			IAssert(keySeen.Len() == nKeys);
		}
	}

public:

	void Report(FILE *f = 0) const
	{
		if (! f) f = stdout;
		const int nBands = hasher.GetNumBands(); 
		fprintf(f, "%s report: %d keys, %d bands\n", typeid(*this).name(), keyHash.Len(), nBands);
		for (int iBand = 0; iBand < nBands; iBand++)
		{
			TStDev stdLens;
			const TBandHash &bandHash = bandHashes[iBand];
			for (int bandCodeId = bandHash.FFirstKeyId(); bandHash.FNextKeyId(bandCodeId); )
			{
				int len = 0; 
				for (int keyId = bandHash[bandCodeId].first; keyId >= 0; 
					keyId = keyHash[keyId].bands[iBand].next) len++;
				stdLens.Add(len);
			}
			fprintf(f, "- Band %d: %d chains; min = %d, max = %d, avg = %.2f +/- %.2f keys\n",
				iBand, stdLens.Count(), (int) stdLens.Min(), (int) stdLens.Max(), stdLens.Avg(), stdLens.Std());
		}
	}

};

//-----------------------------------------------------------------------------
// TBandCodeHasher
//-----------------------------------------------------------------------------
//
// TLshCollection will need to use band codes as keys in hash tables.
// Unfortunately, if we use TB8Set/TB32Set/TBSet as band codes, these classes
// lack the GetPrim/SecHashCd methods.  So here is a template class that
// works around this.  Ideally, GetPrim/SecHashCd should be added directly
// to the bit-set classes.
//
// Note that, technically, TLshCollection uses the THasher::TBandCodeHashFunc
// as the hash function for band codes.  So the TBandCodeHasher is meant to be
// used by concrete THasher implementations to define the TBandCodeHashFunc
// member typedef.

template<typename TBandCode>
struct TBandCodeHasher
{
	static inline int GetPrimHashCd(const TBandCode& Key) { return Key.GetPrimHashCd(); }
	static inline int GetSecHashCd(const TBandCode& Key) { return Key.GetSecHashCd(); }
};

template<>
struct TBandCodeHasher<TB8Set>
{
	static inline int GetPrimHashCd(const TB8Set& Key) { return Key.GetUCh(); }
	static inline int GetSecHashCd(const TB8Set& Key) { return Key.GetUCh(); }
};

template<>
struct TBandCodeHasher<TB32Set>
{
	static inline int GetPrimHashCd(const TB32Set& Key) { int i = (int) Key.GetUInt(); if (i < 0) i = -i; if (i < 0) i = 0; return i; }
	static inline int GetSecHashCd(const TB32Set& Key) { int i = (int) Key.GetUInt(); if (i < 0) i = -i; if (i < 0) i = 0; return i; }
};

// This one is particularly inefficient since TBSet doesn't provide any sort of
// direct access to its underlying vector of ints (TBSet::B4T).
template<>
struct TBandCodeHasher<TBSet>
{
	static inline int GetPrimHashCd(const TBSet& Key) { 
		uint hc = 0, M = 2147483647u; int n = Key.GetBits(); 
		for (int i = 0; i < n; i++) { hc <<= 1; if (Key.GetBit(i)) hc |= 1;
			while (hc >= M) hc -= M; }
		return (int) hc; }
	static inline int GetSecHashCd(const TBSet& Key) { 
		uint hc = 0, M = 1610612741u; int n = Key.GetBits(); 
		for (int i = 0; i < n; i++) { hc <<= 1; if (Key.GetBit(i)) hc |= 1;
			while (hc >= M) hc -= M; }
		return (int) hc; }
};

// TVec's default hash code function simply computes the sum of the
// hash codes of the elements.  In our case this might be problematic if
// the individual hash functions in the band return integers close to 0;
// in that case, the band code is a vector of integers close to 0 and its
// hash code is therefore also an integer still fairly close to 0.
// So here's a hash function that rolls the sum a few bits to the left
// after each addition.
template<typename T>
struct TBandCodeHasher<TVec<T> >
{
	static inline int GetPrimHashCd(const TVec<T>& Key) { 
		uint hc = 0; int n = Key.Len(); 
		enum { LoPart = 13, HiPart = 8 * sizeof(hc) - LoPart };
		for (int i = 0; i < n; i++) 
			hc = ((hc >> LoPart) ^ (hc << HiPart)) + ((uint) Key[i].GetPrimHashCd());  
		return abs(int(hc)); }
	static inline int GetSecHashCd(const TVec<T>& Key) { 
		uint hc = 0; int n = Key.Len(); 
		enum { LoPart = 19, HiPart = 8 * sizeof(hc) - LoPart };
		for (int i = 0; i < n; i++) 
			hc = ((hc >> LoPart) ^ (hc << HiPart)) + ((uint) Key[i].GetSecHashCd());  
		return abs(int(hc)); }
};

//=============================================================================
// Hashers
//=============================================================================

//-----------------------------------------------------------------------------
// IHasher
//-----------------------------------------------------------------------------
//
// This template just shows concisely what a hasher class has to implement in order to be
// used by TLshCollection.  The actual hashers don't have to derive from this
// class or anything of that sort.  
//
// The main role of the hasher is to calculate, for any TDat, a b-tuple
// (where b = GetNumBand()) of "band codes".  Each band code is in turn likely
// to be a r-tuple of hash codes, so that the representation of TDat is ultimately
// based on r*b hash functions; however, the internal composition of a band code
// doesn't matter to the TLshCollection and doesn't need to be exposed by the hasher.
//
// The TLshCollection uses the band codes to perform approximate nearest-neighbor
// queries.  An object x is considered a candidate near neighbor of another object q
// if their band codes are equal in at least one band.  
//
// The band codes are of type TBandCode, which must be exposed by the hasher
// and which must support operator==.  Furtherore, the TLshCollection will be
// using band codes as keys in hash tables, so the hasher has to expose a 
// class named TBandCodeHashFunc, which the TLshCollection will use as a hash function
// in those hash tables.
// 
// The hasher also has to expose three methods: GetNumBands, CalcBandCodes, 
// and BandCodeDistance.
//
// The hasher doesn't actually have to be a template or expose TDat; however,
// the TDat accepted by CalcBandCodes has to match the TDat used by the TLshCollection.
// This allows the same hasher to support multiple TDat types (e.g. by overloading
// CalcBandCodes).

template<typename TDat_>
class IHasher
{
public:
	typedef TDat_ TDat;
	typedef TIntV TBandCode;
	typedef TVec<TBandCode> TBandCodeV;
	typedef TDefaultHashFunc<TBandCode> TBandCodeHashFunc;

	IHasher(TSIn& SIn) { }
	void Save(TSOut &SOut) const { }

	int GetNumBands() const { return -1; }
	void CalcBandCodes(const TDat& dat, TBandCodeV& dest) const;
	// BandCodeDistance should return a measure of how much 'code1' and 'code2' differ, typically from 0 to the
	// number of hash functions in the band.  If code == code2, the return value must be 0,
	// otherwise the return value must be > 0.  This is used by the TLshCollection in cases where
	// the requirement that x and q must match in at least one band code provides too few candidate
	// nearest neighbors; in that case, the search is broadened to include x'es that match one of
	// q's band codes at least approximately (in the sense of having a small BandCodeDistance).
	//
	// This mechanism is based on the assumption that a band code is typically a 
	// r-tuple of hash codes from r different hash functions.  The LshCollection
	// doesn't actually need to know or care about the internal structure of a band code,
	// but the hasher presumably will know about these things and will use them when
	// computing the distance between two band codes.
	int BandCodeDistance(const TBandCode& code1, const TBandCode& code2) const;
};

//-----------------------------------------------------------------------------
// TMinHasher
//-----------------------------------------------------------------------------
//
// Minhashing: the object is a set of integers (of type TItem), represented
// by a vector.  Note that although the order in which items are given in the
// vector doesn't matter as far as the minhasher is concerned, you probably want
// to keep the vectors ordered anyway because the Jaccard distance class (which is 
// likely to be used together with the minhasher) requires the items to be ordered 
// in ascending order.
//
// For a given hash function h, we define minHash(x; h) = min { h(e) : e \in x }.
// Each band i (for i = 1, ..., b) consists of r hash functions h_{i1}, ..., h_{ir}
// and the band code for x is the r-tuple (minHash(x; h_{i1}), ..., minHash(x; h_{ir}).

template<typename TItem_>
class TMinHasher
{
public:
	typedef TItem_ TItem;
	typedef TVec<TItem> TDat;
	enum { BytesPerItem = sizeof(TItem) };
	typedef TIntV TBandCode;
	typedef TVec<TBandCode> TBandCodeV;
	typedef TBandCodeHasher<TBandCode> TBandCodeHashFunc;

protected:
	int nBands, nFuncPerBand;
	TIntV table; 

	// http://en.wikipedia.org/wiki/Tabulation_hashing
	int CalcHashFunc(int iBand, int iFunc, const TItem &item) const
	{
		int result = 0, ofs = (iBand * nFuncPerBand + iFunc) * BytesPerItem * 256;
		for (int iByte = 0; iByte < BytesPerItem; iByte++) {
			int b = int((item >> (8 * iByte)) & 0xff);
			result ^= table[ofs + b];
			ofs += 256; }
		return result;
	}

public:

	TMinHasher(int NBands, int NFuncPerBand) : nBands(NBands), nFuncPerBand(NFuncPerBand) 
	{ 
		TRnd rnd(123);
		table.Gen(nBands * nFuncPerBand * BytesPerItem * 256);
		for (int i = 0; i < table.Len(); i++) {
			int x = rnd.GetUniDevInt(); if (x < 0) { x = -x; if (x < 0) x = 0; }
			table[i] = x; }
	}

	TMinHasher(TSIn& SIn) : nBands(TInt(SIn)), nFuncPerBand(TInt(SIn)), table(SIn) { }
	void Save(TSOut& SOut) const { TInt(nBands).Save(SOut); TInt(nFuncPerBand).Save(SOut); table.Save(SOut); }

	int GetNumBands() const { return nBands; }

	void CalcBandCodes(const TDat& dat, TBandCodeV& dest) const
	{
		if (dest.Len() != nBands) dest.Gen(nBands);
		int nElts = dat.Len();
		for (int iBand = 0; iBand < nBands; iBand++)
		{
			TBandCode &bandCode = dest[iBand];
			if (bandCode.Len() != nFuncPerBand) bandCode.Gen(nFuncPerBand);
			for (int iFunc = 0; iFunc < nFuncPerBand; iFunc++)
			{
				int hMin = 0;
				for (int i = 0; i < nElts; i++) {
					int h = CalcHashFunc(iBand, iFunc, dat[i]);
					if (i == 0 || h < hMin) hMin = h; }
				bandCode[iFunc] = hMin;
			}
		}
	}

	int BandCodeDistance(const TBandCode& code1, const TBandCode& code2) const {
		int d = 0;
		for (int i = 0; i < nFuncPerBand; i++)
			if (code1[i] != code2[i]) d++;
		return d; }
};

//-----------------------------------------------------------------------------
// TRandDirStore
//-----------------------------------------------------------------------------
//
// Stores a set of random directions used by TRandProjHasher.
// Two implementations are provided: vector and hash table.
// See the description of TRandProjHasher for more information.

template<typename TVal_>
class TRandDirStore_Vec
{
public:
	typedef TVal_ TVal;
	typedef TVec<TVal> TValV;
	int nDirs, nDims;
	TValV v; // v[i * nDirs + j] is the i'th component of the j'th direction vector, before normalization
	TValV norms; // norms[j] = sum_i v[i * nDirs + j]^2

protected:

	void InitRow(int dimIdx)
	{
		TRnd rnd(456 + dimIdx);
		for (int i = 0; i < nDirs; i++) {
			double x = rnd.GetNrmDev();
			norms[i] += x * x;
			v[dimIdx * nDirs + i] = x; }
	}

public:

	TRandDirStore_Vec(): nDirs(0), nDims(0) { }
	TRandDirStore_Vec(TSIn& SIn) : nDirs(TInt(SIn)), nDims(TInt(SIn)), v(SIn), norms(SIn) { }
	void Save(TSOut& SOut) const { TInt(nDirs).Save(SOut); TInt(nDims).Save(SOut); v.Save(SOut); norms.Save(SOut); }

	void Init(int nDirs_, int nDims_) {
		nDirs = nDirs_; nDims = nDims_;
		v.Gen(nDims * nDirs); norms.Gen(nDirs); norms.PutAll(0);
		for (int dim = 0; dim < nDims; dim++) InitRow(dim); }

	void EnsureDim(int dimIdx) {
		Assert(dimIdx >= 0);
		if (dimIdx >= nDims) {
			int newDims = dimIdx + 1;
			while (v.Len() < newDims * nDirs) v.Add();
			while (nDims < newDims) InitRow(nDims++); }}

	void AddProj(TValV& dest, int dimIdx, const TVal& val) const {
		Assert(0 <= dimIdx); Assert(dimIdx < nDims);
		for (int i = 0; i < nDirs; i++) {
			double norm = norms[i]; if (norm < 1e-8) norm = 1; else norm = 1.0 / sqrt(norm);
			dest[i] += (val * v[dimIdx * nDirs + i]) * norm; }}

};

template<typename TDimIdx_, typename TVal_>
class TRandDirStore_Hash
{
public:
	typedef TDimIdx_ TDimIdx;
	typedef TVal_ TVal;
	typedef TVec<TVal> TValV;
	typedef THash<TDimIdx, TVoid> TDimIdxH;
	int nDirs, nDims;
	TDimIdxH h; // key: TDimIdx; keyId: 0..nDims-1
	TValV v; // v[i * nDirs + j] is the dims.GetKey(i)'th component of the j'th direction vector, before normalization
	TValV norms; // norms[j] = sum_i v[i * nDirs + j]^2

protected:

	void InitRow(int dimId)
	{
		TRnd rnd(456 + dimId);
		for (int i = 0; i < nDirs; i++) {
			double x = rnd.GetNrmDev();
			norms[i] += x * x;
			v[dimId * nDirs + i] = x; }
	}

public:

	TRandDirStore_Hash() : nDirs(0), nDims(0) { }
	TRandDirStore_Hash(TSIn& SIn) : nDirs(TInt(SIn)), nDims(TInt(SIn)), h(SIn), v(SIn), norms(SIn) { }
	void Save(TSOut& SOut) const { TInt(nDirs).Save(SOut); TInt(nDims).Save(SOut); h.Save(SOut); v.Save(SOut); norms.Save(SOut); }

	void Init(int nDirs_, int nDims_ = 0) {
		nDirs = nDirs_; nDims = 0;
		h.Clr(); v.Clr(); norms.Gen(nDirs); norms.PutAll(0); }

	void EnsureDim(const TDimIdx& dimIdx) {
		int dimId = h.GetKeyId(dimIdx);
		if (dimId >= 0) return; 
		dimId = h.AddKey(dimIdx);
		IAssert(dimId <= nDims);
		if (dimId == nDims) { nDims++; for (int i = 0; i < nDirs; i++) v.Add(); }
		InitRow(dimId); }

	void AddProj(TValV& dest, const TDimIdx &dimIdx, const TVal& val) const {
		Assert(dest.Len() == nDirs);
		const int dimId = h.GetKeyId(dimIdx); 
		Assert(dimId >= 0);
		for (int i = 0; i < nDirs; i++) {
			double norm = norms[i]; if (norm < 1e-8) norm = 1; else norm = 1.0 / sqrt(norm);
			dest[i] += (val * v[dimId * nDirs + i]) * norm; }}

};

//-----------------------------------------------------------------------------
// TBandCodeAccess
//-----------------------------------------------------------------------------
//
// The band code is a r-tuple of hash codes; the following templates
// implement this and provide access to individual components of the tuple.
// Two implementations are provided: vector and bit-set (the latter can
// be used when each compont is +1 or -1).
// See the description of TRandProjHasher for more information.

template<typename TItem_> // the item will typically be a TInt
class TBandCodeAccess_Vec
{
public:
	typedef TItem_ TItem;
	typedef TVec<TItem> TBandCode;

	static void Gen(int nFuncPerBand, TBandCode& dest) { dest.Gen(nFuncPerBand); }
	static void Put(TBandCode& dest, int idx, const TItem item) { dest[idx] = item; }
	static TItem Get(const TBandCode &bandCode, int idx) { return bandCode[idx]; }
};


template<typename TSet_> inline void TBandCodeAccess_Set_Gen(int nFuncPerBand, TSet_& dest) { dest.Clr(); }
template<> inline void TBandCodeAccess_Set_Gen<TBSet>(int nFuncPerBand, TBSet& dest) { dest.Gen(nFuncPerBand); }

template<typename TSet_> // this can be a TB8Set, TB32Set or TBSet -- use the smaller ones if you know that the number of functions per band will be small enough
class TBandCodeAccess_Set
{
public:
	typedef TSet_ TBandCode;

	static void Gen(int nFuncPerBand, TBandCode& dest) { TBandCodeAccess_Set_Gen(nFuncPerBand, dest); }
	static void Put(TBandCode& dest, int idx, int sign) { dest.SetBit(idx, (sign > 0)); }
	static int Get(const TBandCode &bandCode, int idx) { return bandCode.GetBit(idx) ? 1 : -1; }
};

//-----------------------------------------------------------------------------
// TRandProjFinalizer
//-----------------------------------------------------------------------------
//
// Provides postprocessing that convers a projection w^T x into a hash code.
// Two implementations are provided: sign (for use with the cosine distance)
// and floor (for use with the Euclidean distance).
// See the description of TRandProjHasher for more information.

// This converter template can help us avoid some compiler warnings
// that would occur if we cast directly from double to e.g. TInt,
// without an intermediate cast to int.
template<typename TOutVal> struct TRpfConverter { template<typename T> static inline TOutVal Convert(T t) { return (TOutVal) t; } };
template<> struct TRpfConverter<TInt> { template<typename T> static inline TInt Convert(T t) { return (TInt) (int) t; } };
template<> struct TRpfConverter<TSInt> { template<typename T> static inline TSInt Convert(T t) { return (TSInt) (int16) t; } };
template<> struct TRpfConverter<TUInt64> { template<typename T> static inline TUInt64 Convert(T t) { return (TUInt64) (uint64) t; } };
template<> struct TRpfConverter<TUInt> { template<typename T> static inline TUInt Convert(T t) { return (TUInt) (uint) t; } };

template<typename TInVal_, typename TOutVal_>
class TRandProjFinalizer_Floor
{
public:
	typedef TInVal_ TInVal;
	typedef TOutVal_ TOutVal;
	TInVal a;
	TRandProjFinalizer_Floor(TInVal a_) : a(a_) { }
	TRandProjFinalizer_Floor(TSIn& SIn) : a(SIn) { }
	void Save(TSOut& SOut) const { a.Save(SOut); }
	TOutVal Finalize(TInVal projection) const { 
		return TRpfConverter<TOutVal>::Convert(floor(projection / a)); }
};

template<typename TInVal_, typename TOutVal_>
class TRandProjFinalizer_Sign
{
public:
	typedef TInVal_ TInVal;
	typedef TOutVal_ TOutVal;
	TRandProjFinalizer_Sign() { }
	TRandProjFinalizer_Sign(TSIn& SIn) { }
	void Save(TSOut& SOut) const { }
	TOutVal Finalize(TInVal projection) const { return projection > 0 ? 1 : -1; }
};

//-----------------------------------------------------------------------------
// TRandProjHasher
//-----------------------------------------------------------------------------
//
// Here, each object is a (sparse) d-dimensional vector and the hash function h is defined
// as h(x) = sgn(w^T x), where w is a d-dimensional unit vector.
// Each has function has a different w.  If the w's are chosen from a uniform distribution 
// over all unit-length vectors (i.e. the surface of the unit hypersphere) and f = sgn,
// then the resulting family of hash functions will have the property that, for any x and y,
// the probability [over all h] that h(x) != h(y) is proportional to the angle between x and y.
// Thus, the random projection hasher works well with the cosine distance measure.
//
// Choosing unit-length vectors uniformly is a bit tricky.  We can't simply
// start by choosing a vector unifomly from the cube [-1, 1]^d and normalizing it to
// unit length, because this makes "diagonal" vectors more likely than others.
// A better way is to start by choosing d independent numbers from a standardized normal
// distribution, and then normalizing the resulting vector to unit length.  
// 
// If we were interested in the actual values of w^T x instead of just sgn(w^T x), we'd
// have a problem: if we don't know d in advance (e.g. because our vectors are BOW representations
// of documents, so the number of dimensions d can grow as more features are discovered),
// we'd have to occasionally add new components to each w and renormalize it.  This change
// changes w into w' = a w + z, where a is the renormalization factor and z represents the new
// components added to w' to account for the increase in the number of dimensions.
// For an old value of x which has 0 components in these new dimensions, we have z^T x = 0 
// and thus w'^T x = a w^T x + z^T x = a w^T x.  This is different from w^T x, which is the
// has code we computed earlier, when x was first added to the collection (and before the
// number of dimensions has been increased).  However, since the renormalization factor is
// positive, the sign of a w^T x is the same as the sign of w^T x.  Thus, as long as our
// hash codes only keep the sign of each projections, we can easily add new dimensions later
// without invalidating old hash codes.
// - But note that it's important that we do *not* add new dimensions *DURING* the computation
// of the hash codes for an individual x (because this would cause the contributions of
// different source dimensions to be weighted differently, because the norm of the random vector
// changes as additional dimensions are added).  Hence the separation of EnsureDim and AddProj
// in the RandDirStore.
//
// A closely related family of hash functions, h(x) = floor(w^T x  / c) for some real number c > 0,
// can be used together with the Euclidean distance measure.  However, in this case adding new
// dimensions does in fact invalidate old hash codes, so you should provide the correct number of
// dimensions when creating the hasher.
//
// The behavior of the random projection hasher can be customized in several ways 
// using its template parameters:
//
// - TRandDirStore: stores the random directions (the vectors w for all our hash functions).
//   If the dimensions are integers from 0 to d-1 for some reasonable value of d, then it makes 
//   sense to store all the vectors w in a TVec, so use TRandDirStore_Vec.
//   On the other hand, if dimensions are some sort of arbitrary identifiers (e.g. strings; or if there
//   are lots of dimensions that aren't used in any vector that we'll be trying to project),
//   you can use TRandDirStore_Hash which uses a hash table to translate dimension identifiers
//   (TDimIdx) into integers 0..d-1.  In this case, dimension identifiers could even be strings
//   or something else non-numeric.
//
// - TRandProjFinalizer: as we saw above, our hash functions are of the form h(x) = f(w^T x),
//   where f(t) is either sign(t) or floor(t/a).  TRandProjFinalizer implements the function f.
//   You can choose between TRandProjFinalizer_Sign and TRandProjFinalizer_Floor.
//
// - TBandCodeAccess: the band code is a r-tuple of values supplied by our hash functions,
//   and TBandCodeAccess provides access to such a r-tuple.  If you used TRandProjFinalizer_Floor,
//   then each value of each hash function is an integer number, so the r-tuple should be
//   represented by a vector (TBandCodeAccess_Vec).  On the other hand, if you used 
//   TRandProjFinalizer_Sign, then each value of each hash function is just -1 or +1, 
//   so the r-tuple can actually be represented by a bit-set; in this case, use one of
//   TBandCodeAccess_Set<T>, where T = TB8Set, TB32Set or TBSet.  If you have 8 or less
//   hash functions per band, you can use TB8Set; for 32 or less hash functions per band,
//   use TB32Set; finally, TBSet allows an arbitrary number of hash functions per band, 
//   at the cost of higher overheads.  In any case there is little reason to use so many
//   hash functions per band (see the comments at the start of TLshCollection; using too many
//   hash functions per band can lead to very poor performance of the approximate NN queries).

template<typename TDimIdx_, typename TVal_, typename TBandCodeAccess_, typename TRandDirStore_, typename TRandProjFinalizer_>
class TRandProjHasher
{
public:
	// Each of our sparse vectors will be a collection of (TDimIdx, TVal) pairs.
	// Several concrete types of vectors are supported by CalcBandCodes:
	// TVec<TKeyDat<TDimIdx, TVal>>, TVec<TPair<TDimIdx, TVal>>, THash<TDimIdx, TVal>,
	// as well as TVec<TVal>, which is assumed to be a non-sparse vector and the
	// dimension indices are integers from 0 to v.Len()-1.
	typedef TDimIdx_ TDimIdx;
	typedef TVal_ TVal;
	typedef TVec<TVal> TValV;
	typedef TKeyDat<TDimIdx, TVal> TKd;
	typedef TVec<TKd> TKdV;
	typedef TVec<TVal> TValV;
	typedef TPair<TDimIdx, TVal> TPr;
	typedef TVec<TPr> TPrV;
	typedef THash<TDimIdx, TVal> TKdH;

	typedef TBandCodeAccess_ TBandCodeAccess;
	typedef typename TBandCodeAccess::TBandCode TBandCode;

	typedef TVec<TBandCode> TBandCodeV;
	typedef TBandCodeHasher<TBandCode> TBandCodeHashFunc;

	typedef TRandDirStore_ TRandDirStore;
	typedef TRandProjFinalizer_ TRandProjFinalizer;

protected:
	int nBands, nFuncPerBand;
	TRandDirStore randDirStore;
	TRandProjFinalizer finalizer;

public:

	// nDimensions is used to initialize the randDirStore.  If dimensions higher than this
	// later appear in the data, the store will add suitable components to the direction vectors.
	// However, this will lead to wrong results if you use the floor finalizer (and the Euclidean distance)
	// instead of the sign finalizer (and the cosine distance).
	TRandProjHasher(int NBands, int NFuncPerBand, int nDimensions, const TRandProjFinalizer &finalizer_) : nBands(NBands), nFuncPerBand(NFuncPerBand), finalizer(finalizer_)
	{ 
		randDirStore.Init(nBands * nFuncPerBand, nDimensions);
	}

	TRandProjHasher(TSIn& SIn) : nBands(TInt(SIn)), nFuncPerBand(TInt(SIn)), randDirStore(SIn), finalizer(SIn) { }
	void Save(TSOut& SOut) const { TInt(nBands).Save(SOut); TInt(nFuncPerBand).Save(SOut); randDirStore.Save(SOut); finalizer.Save(SOut); }

protected:

	// Given a vector x (in the parameter 'dat'), CalcBandCodes_AddProj computes the 
	// projections w^T x for all nBands * nFuncPerBand hash functions,
	// and stores them in 'projections'.

	void CalcBandCodes_AddProj(const TKdH& dat, TValV &projections) {
		for (int i = dat.FFirstKeyId(); dat.FNextKeyId(i); ) 
			randDirStore.EnsureDim(dat.GetKey(i)); 
		for (int i = dat.FFirstKeyId(); dat.FNextKeyId(i); ) 
			randDirStore.AddProj(projections, dat.GetKey(i), dat[i]); }

	void CalcBandCodes_AddProj(const TKdV& dat, TValV &projections) {
		for (int i = 0; i < dat.Len(); i++) 
			randDirStore.EnsureDim(dat[i].Key); 
		for (int i = 0; i < dat.Len(); i++) 
			randDirStore.AddProj(projections, dat[i].Key, dat[i].Dat); }

	void CalcBandCodes_AddProj(const TPrV& dat, TValV &projections) {
		for (int i = 0; i < dat.Len(); i++) 
			randDirStore.EnsureDim(dat[i].Val1); 
		for (int i = 0; i < dat.Len(); i++) 
			randDirStore.AddProj(projections, dat[i].Val1, dat[i].Val2); }

	// non-sparse vector
	void CalcBandCodes_AddProj(const TValV& dat, TValV &projections) {
		for (int i = 0; i < dat.Len(); i++) 
			randDirStore.EnsureDim(i); 
		for (int i = 0; i < dat.Len(); i++) 
			randDirStore.AddProj(projections, i, dat[i]); }

public:

	template<typename TDat>
	void CalcBandCodes(const TDat& dat, TBandCodeV& dest) 
	{
		if (dest.Len() != nBands) dest.Gen(nBands);
		int nElts = dat.Len();
		// Initialize a vector of projections.
		TValV projections; 
		projections.Gen(nBands * nFuncPerBand);
		projections.PutAll(TVal(0));
		// Calculate the projections w^T x for all hash functions.
		CalcBandCodes_AddProj(dat, projections);
		// Process the projections with the finalizer function and pack them,
		// r at a time, into each of the b band codes.
		for (int iBand = 0; iBand < nBands; iBand++) {
			TBandCodeAccess::Gen(nFuncPerBand, dest[iBand]);
			for (int iFunc = 0; iFunc < nFuncPerBand; iFunc++) {
				TVal proj = projections[iBand * nFuncPerBand + iFunc];
				TRandProjFinalizer::TOutVal projF = finalizer.Finalize(proj);
				TBandCodeAccess::Put(dest[iBand], iFunc, projF); }}
	}

	int GetNumBands() const { return nBands; }

	int BandCodeDistance(const TBandCode& code1, const TBandCode& code2) const {
		int d = 0;
		for (int i = 0; i < nFuncPerBand; i++)
			if (TBandCodeAccess::Get(code1, i) != TBandCodeAccess::Get(code2, i)) d++;
		return d; }
};

//=============================================================================
// Distance functions
//=============================================================================

//-----------------------------------------------------------------------------
// TJaccardDistance
//-----------------------------------------------------------------------------
//
// This distance measure assumes that the objects are sets of TItems,
// and each set is represented by a TVec<TItem> in which the items are
// sorted in ascending order.
//
// The Jaccard distance is defined as J(x, y) = 1 - |x intersection y| / |x union y|.

template<typename TItem_>
class TJaccardDistance
{
public:
	typedef TItem_ TItem;
	typedef TVec<TItem> TDat;
	typedef double TDistance;

	TDistance operator()(const TDat& x, const TDat& y) const {
		int nInter = 0, iy = 0, nx = x.Len(), ny = y.Len();
		for (int ix = 0; ix < nx; ix++) {
			const TItem &xi = x[ix];
			Assert(ix == 0 || x[ix - 1] < xi);
			Assert(iy == 0 || iy >= ny || y[iy - 1] < y[iy]);
			while (iy < ny && y[iy] < xi) iy++;
			if (iy >= ny) break;
			if (xi == y[iy]) { nInter++; iy++; }}
		if (nx + ny <= 0) return 0; else return (nx + ny - 2 * nInter) / double(nx + ny - nInter); }
};

//-----------------------------------------------------------------------------
// TVectorDistanceBase
//-----------------------------------------------------------------------------
//
// A class that can compute dot products between (sparse) vectors consisting of
// (TDimIdx, TVal) pairs.  It will be used as a base class for TCosineDistance
// and TEuclideanDistance.
//
// The vectors can be represented as one of the following:
// - TVec<TVal> - dense vectors (not necessarily of the same length, the missing components of the shorter vector are considered to be 0)
// - TVec<TKeyDat<TDimIdx, TVal> > - sparse vectors (must be sorted in ascending order of TDimIdx)
// - TVec<TPair<TDimIdx, TVal> > - sparse vectors (must be sorted in ascending order of TDimIdx)
// - THash<TDimIdx, TVal> - sparse vectors, represented as hash tables

template<typename TDimIdx_, typename TVal_>
class TVectorDistanceBase 
{
public:
	typedef TDimIdx_ TDimIdx;
	typedef TVal_ TVal;
	typedef TVec<TVal> TValV;
	typedef TKeyDat<TDimIdx, TVal> TKd;
	typedef TVec<TKd> TKdV;
	typedef TPair<TDimIdx, TVal> TPr;
	typedef TVec<TPr> TPrV;
	typedef THash<TDimIdx, TVal> TKdH;
	typedef double TDistance;
public:

	static void Dots(const TValV &x, const TValV &y, double &xx_, double &xy_, double &yy_) 
	{
		double xx = 0, xy = 0, yy = 0;
		int i, nx = x.Len(), ny = y.Len();
		for (i = 0; i < nx && i < ny; i++) {
			double x_ = x[i], y_ = y[i];
			xx += x_ * x_; yy += y_ * y_; xy += x_ * y_; }
		while (i < nx) { double x_ = x[i++]; xx += x_ * x_; }
		while (i < ny) { double y_ = y[i++]; yy += y_ * y_; }
		xx_ = xx; xy_ = xy; yy_ = yy; 
	}

	static void Dots(const TKdV &x, const TKdV &y, double &xx_, double &xy_, double &yy_) 
	{
		double xx = 0, xy = 0, yy = 0;
		int ix = 0, iy = 0, nx = x.Len(), ny = y.Len();
		while (ix < nx || iy < ny) {
			if (ix >= nx) { double y_ = y[iy++].Dat; yy += y_ * y_; }
			else if (iy >= ny) { double x_ = x[ix++].Dat; xx += x_ * x_; }
			else if (x[ix].Key < y[iy].Key) { double x_ = x[ix++].Dat; xx += x_ * x_; }
			else if (y[iy].Key < x[ix].Key) { double y_ = y[iy++].Dat; yy += y_ * y_; }
			else { double x_ = x[ix++].Dat, y_ = y[iy++].Dat; xy += x_ * y_; xx += x_ * x_; yy += y_ * y_; }}
		xx_ = xx; xy_ = xy; yy_ = yy; 
	}

	static void Dots(const TPrV &x, const TPrV &y, double &xx_, double &xy_, double &yy_) 
	{
		double xx = 0, xy = 0, yy = 0;
		int ix = 0, iy = 0, nx = x.Len(), ny = y.Len();
		while (ix < nx || iy < ny) {
			if (ix >= nx) { double y_ = y[iy++].Val2; yy += y_ * y_; }
			else if (iy >= ny) { double x_ = x[ix++].Val2; xx += x_ * x_; }
			else if (x[ix].Val1 < y[iy].Val1) { double x_ = x[ix++].Val2; xx += x_ * x_; }
			else if (y[iy].Val1 < x[ix].Val1) { double y_ = y[iy++].Val2; yy += y_ * y_; }
			else { double x_ = x[ix++].Val2, y_ = y[iy++].Val2; xy += x_ * y_; xx += x_ * x_; yy += y_ * y_; }}
		xx_ = xx; xy_ = xy; yy_ = yy; 
	}

	static void Dots(const TKdH &x, const TKdH &y, double &xx_, double &xy_, double &yy_) 
	{
		double xx = 0, xy = 0, yy = 0;
		for (int ix = x.FFirstKeyId(); x.FNextKeyId(ix); ) {
			double x_ = x[ix]; xx += x_ * x_;
			int iy = y.GetKeyId(x.GetKey(ix));
			if (iy >= 0) xy += x_ * y[iy]; }
		for (int iy = y.FFirstKeyId(); y.FNextKeyId(iy); ) {
			double y_ = y[iy]; yy += y_ * y_; }
		xx_ = xx; xy_ = xy; yy_ = yy; 
	}
}; 

//-----------------------------------------------------------------------------
// TCosineDistance
//-----------------------------------------------------------------------------

template<typename TDimIdx_, typename TVal_>
struct TCosineDistance : public TVectorDistanceBase <TDimIdx_, TVal_>
{
	template<typename T, typename U>
	TDistance operator()(const T& x, const U& y) const {
		double xx, xy, yy; Dots(x, y, xx, xy, yy);
		double norm = xx * yy; 
		// If norm is close to 0, this means x and/or y are 0-length vectors;
		// so in a certain sense, the angle between them (and hence the cosine)
		// is undefined.  But since the purpose of this function is to serve as a
		// distance measure, we'll return 0 if they are both empty, and 1 if
		// just one is nonempty.
		if (norm < 1e-8) return (xx < 1e-8 && yy < 1e-8) ? 0 : 1; 
		norm = 1 / sqrt(norm);
		TDistance cosine = (TDistance) (xy * norm);
		// cos = 1 means identical ---> dist = 0
		// cos = 0 means orthogonal ---> dist = 1
		// cos = -1 means opposite ---> dist = 2
		// Of course, negative cosines will rarely occur since one usually uses 
		// the cosine distance over vectors where all components are >= 0,
		// which ensures that the dot product (and hence the cosine) is also >= 0.
		return 1 - cosine; }
};

//-----------------------------------------------------------------------------
// TEuclideanDistance
//-----------------------------------------------------------------------------

template<typename TDimIdx_, typename TVal_>
struct TEuclideanDistance : public TVectorDistanceBase <TDimIdx_, TVal_>
{
	template<typename T, typename U>
	TDistance operator()(const T& x, const U& y) const {
		double xx, xy, yy; Dots(x, y, xx, xy, yy);
		return (TDistance) (xx + yy - 2 * xy); }
};

//=============================================================================

//-----------------------------------------------------------------------------
// TShingler
//-----------------------------------------------------------------------------
//
// This class converts a string into a set of integers by shingling.
//
// The input string, s, is seen as a sequence of characters: s[0], ..., s[n-1].
// To shingle it, we choose a window length w and compute hash codes of the
// substrings s[0..w-1], s[1..w], s[2..w+1], ..., s[n-w..n-1].
// To compute the hash codes more efficiently, we'll use the Rabin-Karp
// hash function: h(s[i..i+w-1]) = sum_{t=0}^{w-1} b^{w-1-t} s[i+t] mod M,
// where M and b are constant integers.
//
// Several methods are provided:
//
// - ShingleByChars(s, ofs, len, w, sink) is the most generic one.  
//   It assumes that the string we're actually interested in starts at s[ofs]
//   and is 'len' characters long.  It therefore computes shingles for s[i..i+w-1] 
//   for ofs <= i <= ofs + len - w and calls sink(x) for each shingle x.  The only 
//   assumption made about the type of s is that its individual characters can be 
//   accessed as ConvertChar(s[i]).  The call to ConvertChar is used to make
//   sure that individual characters don't end up being interpreted as
//   negative integers.  Thus, s could be TStr, TChA, char *, TUStr, wchar_t *, TIntV, etc.
//
// - ShingleByChars(s, w, dest) is similar to the previous one, but it
//   sets ofs = 0, len = GetLength(s), and the shingles are appended to the 'dest'
//   vector instead of being passed to a sink function.  GetLength() is an overloaded
//   method function which either calls s.Len() or strlen/wcslen, depending on the type of s.
// 
// - ShingleString(s, ofs, len) computes just one shingle, namely for the
//   entire string s[ofs...ofs+len-1], and returns it as the return value of the function.
//
// - ShingleString(s) is like the previous one, with ofs = 0 and len = GetLength(s).
//
// - ShingleByWords(words, ofs, len, w, sink) computes one shingle (using ShingleString)
//   for each word from words[ofs] to words[ofs + len - 1].  The resulting vector of
//   'len' shingles is then itself shingled (using ShingleByChars) with a window length 'w'.
//   Thus, the sink receives one shingle for each w-gram of words from the original document.
//   The type of 'words' can be anything that supports operator [] and returns strings
//   of a type that ShingleString can work with.
//
// - ShingleByWords(words, dest) is like the previous one, with ofs = 0 and len = words.Len().
//
// Note that no functionality to break a string into words is provided. 
// This should be done by the user, who is in a better position to determine what exactly
// counts as a word, whether to remove stopwords, etc.  Word breaking is already available 
// elsewhere in glib, e.g. TStr::SplitOnWs, TStr::SplitOnNonAlNum, TUStr::GetWordUStrV.
//
// Shingling a string gives you a list of integers, which you can use as a TDat
// in a TLshCollection (together with the MinHasher and the Jaccard distance measure),
// but note that you should sort each such list first before it can be used with the 
// Jaccard measure.

template<
	typename TShingle_ = int,
	typename TInternalInt_ = long long
>
class TShingler
{
public:
	typedef TShingle_ TShingle;
	typedef TInternalInt_ TInternalInt;
	typedef TVec<TInternalInt> TInternalIntV;
protected:
	TInternalInt M, b;
public:
	// TShingle must be able to contain any of the values 0..M-1.
	// TInternalInt must be able to contain any of the values 0..(M * max(2, b) - 1).
	TShingler() : M((TInternalInt) 2147483647), b((TInternalInt) 17) { }
	TShingler(const TInternalInt M_, const TInternalInt b_) : M(M_), b(b_) { }

	template <typename TVector>
	class TVectorSink
	{
	protected:
		TVector &v;
	public:
		TVectorSink(TVector &V) : v(V) { }
		template <typename TItem> void operator()(const TItem &item) { v.Add(item); }
	};

	template <typename T> static inline TInternalInt ConvertChar(T t) { return t; }
	template <> static inline TInternalInt ConvertChar(char c) { return (uchar) c; }
	template <> static inline TInternalInt ConvertChar(short c) { return (ushort) c; }
	template <> static inline TInternalInt ConvertChar(TSInt c) { return (ushort) c.Val; }

	template <typename TString, typename TOffset, typename TSink>
	void ShingleByChars(const TString &s, TOffset ofs, TOffset len, TOffset window, TSink &sink)  const
	{
		if (window > len) return;
		// hc(s[i..i+w-1]) = sum_{t=0}^{w-1} b^{w-1-t} s[i+t]    mod M.
		// Therefore, hc(s[i+1..i+w]) = (hc(s[i..i+w-1]) - b^{w-1} s[0]) * b + s[i+w].
		TInternalInt cur = 0, bw1 = 1;
		for (TOffset i = 0; i < window - 1; i++) {
			bw1 = (bw1 * b) % M; 
			cur = (((cur * b) % M) + (ConvertChar(s[ofs + i]) % M)) % M; }
		if (bw1 != 0) bw1 = M - bw1;
		// At this point, bw1 = -b^{w-1} mod M, and cur = sum_{t=0}^{w-2} b^{w-1-t} s[t] mod M.
		for (TOffset i = ofs + window - 1; i < ofs + len; i++) {
			// At this point, cur = sum_{t=0}^{w-2} b^{w-1-t} s[i - w+1 + t] mod M.
			cur = ((cur * b) % M + ConvertChar(s[i])) % M;
			// At this point, cur = sum_{t=0}^{w-1} b^{w-1-t} s[i - w+1 + t] mod M.
			sink((TShingle) cur);
			// Remove the leftmost term from 'cur'.
			cur = (cur + (bw1 * (ConvertChar(s[i - window + 1]) % M)) % M) % M; }
	}

	// An simpler and less efficient version of ShingleByChars.  It's intended for 
	// testing only; it should always produce exactly the same results as ShingleByChars.
	template <typename TString, typename TOffset, typename TSink>
	void ShingleByChars_Slow(const TString &s, TOffset ofs, TOffset len, TOffset window, TSink &sink)  const
	{
		len += ofs;
		while (ofs + window <= len)
		{
			TInternalInt cur = 0;
			for (TOffset i = 0; i < window; i++) {
				cur = (cur * b) % M;
				cur = (cur + (ConvertChar(s[ofs + i]) % M)) % M; }
			sink((TShingle) cur); 
			ofs++;
		}
	}

public:

	static inline int GetLength(const char *s) { return s ? (int) strlen(s) : 0; }
	static inline int GetLength(const uchar *s) { return s ? (int) strlen((const char *) s) : 0; }
	static inline int GetLength(const wchar_t *s) { return s ? (int) wcslen(s) : 0; }
	static inline int GetLength(const TStr& s) { return s.Len(); }
	static inline int GetLength(const TChA& s) { return s.Len(); }
	static inline int GetLength(const TUStr& s) { return s.Len(); }
	static inline int GetLength(const short *s) { int l = 0; while (s[l]) l++; return l; }
	template<typename T> static inline int GetLength(const TVec<T>& s) { return s.Len(); }

	template <typename TString, typename TVector>
	void ShingleByChars(const TString& s, int window, TVector& dest, bool clrDest = true) const {
		if (clrDest) dest.Clr(); ShingleByChars(s, 0, GetLength(s), window, TVectorSink<TVector>(dest)); }

	template <typename TString, typename TOffset>
	TShingle ShingleString(const TString &s, TOffset ofs, TOffset len) const {
		TInternalInt cur = 0; for (TOffset i = 0; i < len; i++) 
			cur = (((cur * b) % M) + ConvertChar(s[ofs + i])) % M; 
		return (TShingle) cur; }

	template <typename TString>
	TShingle ShingleString(const TString &s) const { return ShingleString<TString, int>(s, 0, GetLength(s)); }

	template <typename TVector, typename TOffset, typename TSink>
	void ShingleByWords(const TVector &words, TOffset ofs, TOffset len, TOffset window, TSink &sink) const
	{
		if (ofs + window > len) return;
		TInternalIntV wordShingles; wordShingles.Gen(len - ofs);
		for (TOffset i = ofs; i < len; i++) {
			TInternalInt cur = (TInternalInt) ShingleString(words[i]);
			cur = (cur * b) % M; // imagine a null separator between words, so that a shingle of "ab c" is different from a shingle of "a bc"
			wordShingles[i - ofs] = cur; }
		ShingleByChars(wordShingles, 0, len - ofs, window, sink);
	}

	template <typename TWordVector, typename TShingleVector, typename TOffset>
	void ShingleByWords(const TWordVector &words, TOffset window, TShingleVector &dest, bool clrDest = true) const {
		if (clrDest) dest.Clr(); ShingleByWords(words, 0, words.Len(), window, TVectorSink(dest)); }
};

//-----------------------------------------------------------------------------
// Shortcuts for easier declaration/construction of LSH classes
//-----------------------------------------------------------------------------
//
// The following templates provide easier ways to declare and construct
// TLshCollection instances with typical combinations of template parameters.
// Each of them defines TDat, THasher, TDistFunc, TLsh and PLsh.
// A convenient construction function (New) is also provided.
// The key type (TKey) is left as a template parameter.
//
// Example:
//    Lsh_IntSet_Jaccard<TStr>::PLsh lsh = Lsh_IntSet_Jaccard::New(nBands, nFuncPerBand);


template<typename TKey>
struct Lsh_IntSet_Jaccard 
{
	typedef TInt TItem;
	typedef TVec<TItem> TDat;
	typedef TMinHasher<TItem> THasher;
	typedef TJaccardDistance<TItem> TDistFunc;
	typedef TLshCollection<TKey, TDat, THasher, TDistFunc> TLsh;
	typedef TPt<TLsh> PLsh;
	static PLsh New(int nBands, int nFuncPerBand) { return new TLsh(THasher(nBands, nFuncPerBand)); }
};

template<typename TKey>
struct Lsh_IntFltKdV_Cosine_8 // max 8 functions per band
{
	typedef TInt TDimIdx;
	typedef TFlt TVal;
	typedef TVec<TKeyDat<TDimIdx, TVal> > TDat;
	typedef TRandProjFinalizer_Sign<TVal, int> TFinalizer;
	typedef TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Set<TB8Set>, TRandDirStore_Vec<TVal>, TFinalizer> THasher;
	typedef TLshCollection<TKey, TDat, THasher, TCosineDistance<TInt, TFlt> > TLsh;
	typedef TPt<TLsh> PLsh;
	static PLsh New(int nBands, int nFuncPerBand) { IAssert(nFuncPerBand <= 8); 
		return new TLsh(THasher(nBands, nFuncPerBand, 0, TFinalizer())); }
};

template<typename TKey>
struct Lsh_IntFltKdV_Cosine_32 // max 32 functions per band [NOTE: it's very unlikely that you actually want that many, let alone more!]
{
	typedef TInt TDimIdx;
	typedef TFlt TVal;
	typedef TVec<TKeyDat<TDimIdx, TVal> > TDat;
	typedef TRandProjFinalizer_Sign<TVal, int> TFinalizer;
	typedef TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Set<TB32Set>, TRandDirStore_Vec<TVal>, TFinalizer> THasher;
	typedef TLshCollection<TKey, TDat, THasher, TCosineDistance<TInt, TFlt> > TLsh;
	typedef TPt<TLsh> PLsh;
	static PLsh New(int nBands, int nFuncPerBand) { IAssert(nFuncPerBand <= 32); 
		return new TLsh(THasher(nBands, nFuncPerBand, 0, TFinalizer())); }
};

template<typename TKey>
struct Lsh_IntFltKdV_Cosine_X // arbitrary number of functions per band
{
	typedef TInt TDimIdx;
	typedef TFlt TVal;
	typedef TVec<TKeyDat<TDimIdx, TVal> > TDat;
	typedef TRandProjFinalizer_Sign<TVal, int> TFinalizer;
	typedef TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Set<TBSet>, TRandDirStore_Vec<TVal>, TFinalizer> THasher;
	typedef TCosineDistance<TDimIdx, TVal> TDistFunc;
	typedef TLshCollection<TKey, TDat, THasher, TDistFunc> TLsh;
	typedef TPt<TLsh> PLsh;
	static PLsh New(int nBands, int nFuncPerBand) { 
		return new TLsh(THasher(nBands, nFuncPerBand, 0, TFinalizer())); }
};

template<typename TKey>
struct Lsh_IntFltKdV_Euclidean
{
	typedef TInt TDimIdx;
	typedef TFlt TVal;
	typedef TVec<TKeyDat<TDimIdx, TVal> > TDat;
	typedef TInt TItem;
	typedef TRandProjFinalizer_Floor<TVal, TItem> TFinalizer;
	typedef TRandProjHasher<TDimIdx, TVal, TBandCodeAccess_Vec<TItem>, TRandDirStore_Vec<TVal>, TFinalizer> THasher;
	typedef TEuclideanDistance<TDimIdx, TVal> TDistFunc;
	typedef TLshCollection<TKey, TDat, THasher, TDistFunc> TLsh;
	typedef TPt<TLsh> PLsh;
	// nDimensions = the dimensionality of the space from which our feature vectors will be taken.
	// projectionGranularity = each band code is a tuple of hash codes of the form h(x) = floor(w^T x / projectionGranularity),
	//   for various random directions x.
	static PLsh New(int nBands, int nFuncPerBand, int nDimensions, double projectionGranularity) { 
		return new TLsh(THasher(nBands, nFuncPerBand, nDimensions, TFinalizer(projectionGranularity))); }
};

#endif // ____LSH_H_INCLUDED____