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


/*.................................................................................................................................................................*/

// function to check a given number n is prime or not:
// T.C. :- O(sqrt(n))
bool is_prime(int n){
    if(n == 0 || n == 1) return 0;
    if(n == 2 || n == 3) return 1;
    if(n%2 == 0 || n%3 == 0) return 0;
    for(int i = 5; i*i <= n; i += 6){
        if(n%i == 0 || n%(i+2) == 0) return 0;
    }
    return 1;
}

/*.................................................................................................................................................................*/

// we can atmax create an array/vector of size 10^6 inside our main function.
// when we create an array/vector globally everything is initialised to 0 by default.
// to create an array/vector of size 10^7 of type int we need to declare it globally i.e outside main function.
// we can create an array/vector of size 10^8 of type bool globally.

// sieve of eratosthenes -> given number is prime or not for multiple queries.
// T.C. :- O(n log logn)
void sieve_of_eratosthenes(vector<int> &is_prime, int n){
    is_prime[0] = is_prime[1] = 0;
    for(int i = 2; i*i <= n; i++){
        if(is_prime[i]){
            for(int j = i*i; j <= n; j+=i){
                is_prime[j] = 0;
            }
        }
    }
}

int main(){
    int n = 1e6;
    vector<int> is_prime(n+1, 1);

    sieve_of_eratosthenes(is_prime, n);

    // loop for printing output
    for(int i = 0; i <= n; i++){
        cout << i << " " << is_prime[i] << endl;
    }

    return 0;
}

                                            ------------------------------------------------------------------------

// sieve of eratosthenes -> stores the smallest prime factor for a given number in sieve vector.
// T.C. :- O(n log logn)						    
void sieve_of_eratosthenes(vector<int> &sieve, int n){
    for(int i = 2; i*i <= n; i++){
        if(sieve[i] == i){
            for(int j = i*i; j <= n; j+=i){
                if(sieve[j] == j) sieve[j] = i;
            }
        }
    }
}

int main(){
    int n = 1e6;
    vector<int> sieve(n+1, 1);
    for(int i = 0; i <= n; i++) sieve[i] = i; // store number from 1 to n in vector

    sieve_of_eratosthenes(sieve, n);

    // prints smallest prime factor for every number
    for(int i = 0; i <= n; i++){
        cout << i << " " << sieve[i] << endl;
    }

    return 0;
}

/*.................................................................................................................................................................*/

// function to print prime factorization of a number n:

// Method-1:
// T.C. :- O(sqrt(n)) - good if we have just one query.                           
void prime_factorization(int n){
    for(int i = 2; i*i <= n; i++){
        while(n%i == 0){
            cout << i << " ";
            n = n/i;
        }
    }
    if(n >= 2) cout << n << " ";
}

int main(){
    int n;
    cin >> n;
    prime_factorization(n);

    return 0;
}

                                            ------------------------------------------------------------------------

// Method-2:
// T.C. :- O(n log logn) - more suited if we have a lot of queries as we just create sieve once and then just O(logn) time is required to print prime factorization                           
void sieve_of_eratosthenes(vector<int> &sieve, int n){
    for(int i = 2; i*i <= n; i++){
        if(sieve[i] == i){
            for(int j = i*i; j <= n; j+=i){
                if(sieve[j] == j) sieve[j] = i;
            }
        }
    }
}

void prime_factorization(int n){
    vector<int> sieve(n+1, 1);
    for(int i = 0; i <= n; i++) sieve[i] = i; // store number from 1 to n in vector

    sieve_of_eratosthenes(sieve, n);

    // printing prime factorization:
    while(n != 1){
        cout << sieve[n] << " ";
        n = n/sieve[n];
    }

    return;
}

int main(){
    int n;
    cin >> n;
    
    prime_factorization(n);

    return 0;
}

/*.................................................................................................................................................................*/
















