// IndexStructures.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

namespace BtreeTest 
{
	int BtreeTestMain(int argc, char** argv);
}

namespace LshTest 
{
	int LshTest();
	int LshTest2();
	int LshTest3();
	int ShingleTest();
	int ShingleTest2();
	void LshExample();
	void LshExample2();
}


int main(int argc, char** argv)
{
	TStr s = "Hello World!";
	printf("%s\n", s.CStr());
	if (false) return BtreeTest::BtreeTestMain(argc, argv);
	if (false) LshTest::LshTest();
	if (true) LshTest::LshTest2();
	if (false) LshTest::LshTest3();
	if (false) LshTest::ShingleTest();
	if (false) LshTest::ShingleTest2();
	if (false) { LshTest::LshExample(); LshTest::LshExample2(); }
	return 0;
}
