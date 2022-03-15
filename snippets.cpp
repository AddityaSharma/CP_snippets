#include<bits/stdc++.h>
using namespace std;

int mod = 1e9+7;

// modular exponentiation: -> computing (base)^n
long long power(int base, int n){
	long long ans = 1;
	while(n != 0){
		if(n%2) {
			n = n-1;
			ans = (ans * base) % mod;
		}else{
			n = n/2;
			base = (base * base) % mod;
		}
	}
	return ans;
}


// modulo division: -> computing c = (b/a) % mod 
long long power(int base, int n){
	long long ans = 1;
	while(n != 0){
		if(n%2) {
			n = n-1;
			ans = (ans * base) % mod;
		}else{
			n = n/2;
			base = (base * base) % mod;
		}
	}
	return ans;
}

int main(){
	int mod = 1e9+7;
	int a = 6;
	int b = 7;
	int c = (b * power(a, mod-2)) % mod; // c = (b/a)%mod.
	cout << c << endl;
	return 0;
}
