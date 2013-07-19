#ifndef ____TESTCOMMON_H____
#define ____TESTCOMMON_H____

#define TTimer_DEFINED
class TTimer
{
public:
	unsigned long long total, since;
	int running;

	inline static unsigned long long Now() {
		FILETIME ct, et, kt, ut;
		GetProcessTimes(GetCurrentProcess(), &ct, &et, &kt, &ut);
		ULARGE_INTEGER ul; ul.LowPart = ut.dwLowDateTime;
		ul.HighPart = ut.dwHighDateTime;
		return ul.QuadPart; }

	TTimer() { Clr(); }
	void Clr() { total = 0; since = 0; running = 0; }
	void Start() { if (++running == 1) since = Now(); }
	TTimer &Stop() { if (--running == 0) total += Now() - since; return *this; }
	unsigned long long Ticks() const {
		unsigned long long t = total;
		if (running > 0) t += Now() - since;
		return t; }
	double Sec() const { return double(Ticks()) / 10000000.0; }
	TTimer &Add(unsigned long long ticks) { total += ticks; return *this; }
	TTimer &Add(const TTimer& other) { total += other.Ticks(); return *this; }
};

#endif // ____TESTCOMMON_H____