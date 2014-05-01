#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <algorithm>
#include <sstream>
#include <ctime>
#include <bitset>

typedef long long ll;

using namespace std;

bool printOneNonePrime;
bool printOnePrime;
bool printRanTrace;
bool printPowModTrace;
int lineNum;
// for unkown bit length numbits = -1;
string printBin(int n,int numbits)
{
	int c = n;
	int t;
	string bitstr = "";
	if(numbits > 0){ // if bit length is known
		for(int i=numbits;i>0;i--){
			t = 1;
			t <<= i-1; 
			if(c&t){
				bitstr += "1";
			}else{
				bitstr += "0";
			}
		}
	}else{ // if bit length is unknown;
		while(c>0){
			if(c&1){
				bitstr = "1" + bitstr;
			}else{
				bitstr = "0" + bitstr;
			}
			c >>= 1;
		}
	}
	return bitstr;
}

string vecToStr(vector<bool> v)
{
	string ret= "";
	for(int i=0;i<v.size();i++)
	{
		ret += ((v[i]) ? "1" : "0");
	}
	return ret;
}

void printKeyMap(map<string,int> m,int primeBitLen){
	cout << "p is " << m["p"] << "; bits for p " << printBin(m["p"],primeBitLen) << endl;
	cout << "q is " << m["q"] << "; bits for q " << printBin(m["q"],primeBitLen) << endl;
	cout << "n is " << m["n"] << "; bits for n " << printBin(m["n"],-1) << endl;
	cout << "e is " << m["e"] << "; bits for e " << printBin(m["e"],-1) << endl;
	cout << "d is " << m["d"] << "; bits for d " << printBin(m["d"],-1) << endl; 
}

// multiply modulo n
ll mulmod(ll a,ll b,ll n)
{
	ll x = 0;
	ll y = a % n;
	while(b){
		if(b&1)
		{
			x = (x + y) % n;
		}
		y <<= 1; // y multiply by 2;
		b >>= 1; // b divided by 2;
	}
	return x;
}

// fast exponentiation
ll powmod(ll a,ll x,ll n)
{
	if(printPowModTrace){
		cout << "powmod trace "<< a << " to the power of "<< x << " modulo " << n << endl;
	}

	ll y = 1;
	int k = floor(log2(x));
	for(int i=k;i>=0;i--)
	{
		int t = 1;
		t <<= i;
		y = mulmod(y,y,n);
		if(x&t){ // if the ith bit is 1;
			y = mulmod(a,y,n);
		}
		if(printPowModTrace){
			cout << y << endl;
		}
	}
	return y;
}

//generate n number of random bits; v ==> verbose
int ranBits(int n,bool v)
{
	int ret = 0;
	srand(time(NULL));

	if (!printRanTrace) 
	{
		cout << "\n#####line 97######\n" << endl;
	}

	for(int i=0;i<n;i++)
	{
		int ranNum = rand();
		if (!printRanTrace) 
		{
			cout << "random number is " << ranNum << " ";
		}

		if( ranNum & 1 )
		{
			if(!printRanTrace)
			{
				cout << " random bit is 1 "<<endl;  
			}
			ret += ((2 << i) >> 1);
		}else{
			if(!printRanTrace)
			{
				cout << " random bit is 0 "<<endl; 
			} 
		}
	}
	if(!printRanTrace) 
	{
		cout << n <<" random bits in binary form is "<< printBin(ret,n)<<endl; 
	}
	printRanTrace = true;
	return ret;
}

//generate random Number with n bits; v == verbose option
int getRanNum(int n,bool v)
{
	int middle = ranBits(n-2,v);
	middle <<= 1;
	middle++; //add the least significant bit;
	int mostSig = 1;
	mostSig <<= n-1;
	middle += mostSig; // add the most significant bit;
	return middle;
}

// Miller Rabin for a single a where 0< a <= n - 1 
bool millerRabin(int a,int x,int n)
{
	int k = (int)floor(log2((double)x));
	int y = 1;
	int z ;
	for(int i=k;i>=0;i--)
	{
		z = y;
		y = mulmod(y,y,n);
		int t = 1;
		t <<= i;
		if( y == 1 && z != 1 && z != n-1 )
		{
			return false;
		}

		if(x&t)
		{
			y = mulmod(y,a,n);
		}
	}
	if(y!=1){
		return false; 
	}else{
		return true; // a is perhaps prime
	}
	return 0;
}

// test number p with x random numbers a (0 <= a <= p-1); v => verbose
bool primalityTest(int p, int x){
	srand(time(NULL));
	int a;
	for(int i=0;i<x;i++){
		a = rand() % p;
		a = (a==0) ? ++a : a; // a cannot be zero
		if(!millerRabin(a,p-1,p))
		{	
			if(!printOneNonePrime){
				cout << "\n#####line 112######\n" << endl;
				cout << "n " << p << " is not a prime; number a used to test is " << a << endl;
				printOneNonePrime = true;
			}
			return false;
		}
	}

	if(!printOnePrime){
		cout << "\n#####line 115######\n" << endl;
		cout << "n " << p << " is probably a prime; one number a used to test is " << a << endl;
		printOnePrime = true;
	}

	return true;
}

// extended euclidean algorithm
int getInverse(int a, int b)
{
	int q = a/b;
	int r = a%b;
	int s[3];
	int t[3];
	int j = 2;
	int n = b;
	s[0] = 0;
	t[0] = 1;
	s[1] = 1;
	t[1] = q*-1;

	while(b%r != 0)
	{
		a = b;
		b = r;
		q = a /b;
		r = a % b;
		s[j%3] = s[(j-2)%3] + (q * -1 * s[(j-1)%3]);
		t[j%3] = t[(j-2)%3] + (q * -1 * t[(j-1)%3]);
		j++;
	}

	int ret = s[(j-1)%3];
	return ((ret < 0 ) ? ret + n : ret);
}

//find greatest common divisor using euclidean algorithm
int gcd(int a,int b)
{
	if(a == b)
		return b;

	if(a==0 || b == 0)
		return -1;

	int i = 1;
	int q = a/b;
	int r = a%b;	

	while(r!=0){
		a = b;
		b = r;
		q = a/b;
		r = a%b;
		i++;
	}

	return b;
}

map<string,int> RSA(int primeBitLen)
{
	// number of bits for prime number;
    int p = -1;
	int q = -1;
	int temp,e,d,phi_n,n;
	int c=0;
	bool v = false;

	do{
		do{
			for(int i=0;i<2;i++){
				do{	
				v = (c==0) ? true : false;
				temp = getRanNum(primeBitLen,v);   // sverbose only when c == 0
				c++;
				}while(!primalityTest(temp,20));
				switch(i){
					case 0:
						p = temp;
						break;
					case 1:
						q = temp;
						break;
				}
			}
		}while(q == p); 

		cout << "q is " << q << " and p is " << p << endl;
		cout << "q " << printBin(q,primeBitLen)<< endl;
		cout << "p " << printBin(p,primeBitLen) << endl;

		cout << "\n#####line 133######\n" << endl;
		e = 3;
		phi_n = (p-1)*(q-1);
		n = p*q;
		while(gcd(phi_n,e)!=1){
			cout << "trying e with value " << e << endl;
			e++;
		}
	}while(e >= phi_n);
	d = getInverse(e,phi_n);
	cout << "\n#####line 143######\n" << endl;
	cout << "d is " << d << endl;

	map<string,int> m;
	m["p"] = p;
	m["q"] = q;
	m["n"] = n;
	m["e"] = e;
	m["d"] = d;

	return m;
}

/* 
* generate a random number u from n such that 
* the index position of the most significant bit in n is k
* and the index position of the most significant bit in u is k-1;
*/
int getRanNumFromN(int n)
{
	int k = floor(log2(n));
	int u = 1;
	u <<= k-1;
	int remain = getRanNum(k-1,false);
	u += remain;
	cout << "\n#####line 195######\n" << endl;
	cout << "k " << k << " u " << u << endl;

	cout << "\n#####line 197######\n" << endl;
	cout << "u as bits " << printBin(u,32)<< endl;
	return u;
}

vector<bool> genR(string s,int n,int e)
{
	vector<bool> ret (14*8,0);
	const char * c = s.c_str();
	// first 6 bytes
	int diff = 6-s.length(); // account for the blank
	int a;
	for(int i=0;i<6;i++)
	{
		if(i<s.length()){
			a = c[i];
			int bitIndex = 0;
			while(a!=0){	
				ret[i*8 - 1 + 8 - bitIndex] = ((a&1)?1:0);
				a >>= 1;
				bitIndex++;
			}
		}
	}

	// bytes 7 to 10
	a = n;
	for(int i=6*8;i<10*8;i++)
	{
		ret[80 + 6*8 -1-i] = ((a&1)? 1 : 0);
		a >>= 1;
	}

	// bytes 11 to 14
	a = e;
	for(int i=10*8;i<14*8;i++)
	{
		ret[14*8 + 10*8 -1-i] = ((a&1)? 1 : 0);
		a >>= 1;
	}

	return ret;
}

// hash alice's pair r to generate a digital certificate
vector<bool> hashR(vector<bool> r)
{
	vector<bool> ret (8,-1);

	if(r.size() != 14*8){
		return ret;
	}

	for(int j=0;j<8;j++){
		ret[j] = r[j];
		for(int i=1;i<14;i++){
			ret[j] = ret[j]^r[j+i*8];
		}
	}
	return ret;
}

// hash u sent from Bob or v sent from Alice
bitset<8> hashUV(int u)
{
	bitset<32> b1 (u);
	bitset<8> b2;
	for(int j=0;j<8;j++){
		b2[j] = b1[j];
		for(int i=1;i<4;i++)
		{
			b2[j] = b2[j] ^ b1[j+i*8];
		}
	}
	return b2;
}


// converts a vector<bool> into an int
int vectToInt(vector<bool> h)
{
	int ret = 0;
	for(int i=0;i<8;i++)
	{
		if(h[i])
			ret += (int)pow((double)2,(double)h.size()-1-i);
	}
	return ret;
}

int main()
{
	lineNum = 0;

	// get p and q;
	cout << "\n\n=============Question 4 Building a RSA System==============" << endl;
	int primeBitLen = 7;
	map<string,int> m = RSA(primeBitLen);
	cout << "Found e !" << endl;
	cout << "Keys generation complete !"<< endl;
	cout << "\n#####line 147######\n" << endl;
	printKeyMap(m,primeBitLen);
	cout << "\n\n=============Question 5 Creating a Digital Certificate==============" << endl;
	cout << "Generating keys for Trent " << endl;
	map<string,int> mt = RSA(primeBitLen);
	printKeyMap(mt,primeBitLen);
	cout << "Generating keys for Alice " << endl;
	map<string,int> ma = RSA(primeBitLen);
	printKeyMap(ma,primeBitLen); 
	vector<bool> r = genR("Alice",ma["n"],ma["e"]);
	vector<bool> h_v = hashR(r);
	int h_int = vectToInt(h_v);
	int s = powmod(h_int,mt["d"],mt["n"]);
	cout << "\n#####line 175######\n" << endl;
	cout << "r " << vecToStr(r) << " h(r) " << vecToStr(h_v) << " s " << printBin(s,32) << endl; 
	cout << "\n#####line 177######\n" << endl;
	cout << "h(r) as integer " << h_int << " s as integer " << s << endl;
	cout << "\n\n=============Question 6 Alice authenticates herself to Bob ==============" << endl;

	// Bob picks a random number u
	int u = getRanNumFromN(ma["n"]);

	
	// Alice computes hu
	bitset<8> h_u = hashUV(u);
	
	// Alice encrytp it with her private key
	int v =  powmod(h_u.to_ulong(),ma["d"],ma["n"]);

	printPowModTrace = true;
	//calculation of E(e,v)
	cout << "\n#####line 207######\n" << endl;
	int e_v = powmod(v,ma["e"],ma["n"]);

	cout << "\n#####line 204######\n" << endl;
	cout << "u, h(u),v,E(e,v) as bits: "<< endl;
	cout << "u " << printBin(u,32) <<" h_u " <<h_u.to_string() << " v " << printBin(v,32) << " E(e,v) " << printBin(e_v,32) << endl;

	cout << "u, h(u),v,E(e,v) as integers: "<< endl;
	cout << "u " << u <<" h(u) " <<h_u.to_ulong() << " v " << v << " E(e,v) " << e_v << endl;

	/*
	for(int i=0;i<h_v.size();i++)
	{
		cout << h_v[i];
	}
	cout << endl;
	
	int h_u = powmod(u,ma["d"],ma["n"]);  */

	return 0;
}