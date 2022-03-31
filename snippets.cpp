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

// segmented sieve:
// in range [l,r] -> tell number of primes in that range.
// constrains : 1) 1 <= l <= r <= 1e10
//              2) 1 <= r-l <= 1e3

// step-1: store all the primes that will be responsible for marking primes in l to r, i.e all primes till sqrt(r).
// step-2: make array of size (r-l), where 0th index will be representing (0+l)th index   
// step-3: mark all multiples of primes stored in step-1 in the range [l,r] 
// step-4: count elements marked with true.

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

void segmented_sieve(){
    int n = 1e5;
    vector<int> is_prime(n, 1);
    sieve_of_eratosthenes(is_prime, n);

    int l, r; // range - [l..r]
    cin >> l >> r;

    // collect all primes till sqrt(r).
    vector<int> primes;
    for(int i = 2; i*i <= r; i++){
        primes.push_back(i);
    }

    // create dummy array of size (r-l)
    vector<int> dummy((r-l)+1, 1);

    // mark all multiples of primes:
    for(auto pr : primes){
        int first_multiple = (l / pr) * pr;
        if(l % pr) first_multiple += pr; 
        for(int j = max(first_multiple, pr*pr); j <= r; j += pr){
            dummy[j-l] = 0;
        }
    }

    // iterate and count primes:
    int count = 0;
    for(auto x : dummy){
        if(x == 1) count++;
    }

    cout << count;
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

// Bit Manipulation Tricks: 
// check for odd : if(n & 1) n is odd. i.e (n % 2) == 1 is equivalent to n&1 == 1
// check for even : if(!(n & 1)) n is even. (n % 2) == 0 is equivalent to n&1 == 0
// right shift operator : >> - divide by 2
// left shift operator : << - multiply by 2
// pow(2,x) -> 1 << x

// int : 32 bits
// long long : 64 bits
// most significant bit(MSB) : log2(n);

// iterating over bits of a given number 'n' :
while(n){
    int bit = n & 1; // bits are traversed from r to l
    n = n >> 1;
}

// checking ith bit is set or not : if(n & (1 << i) != 0) -> ith bit is set
// set a bit at ith position : n = n | (1<<i)
// count of set bits : __builtin_popcount(n); n -> interger
// count of set bits : __builtin_popcountll(n); n -> long long

// xor : even number of 1's is 0 -> imp way to find xor.
// Q : for given n, find xor from 1 to n, n <= 1e18 -> imp ques.
// observations : n % 4 == 0 -> xor = n
//                n % 4 == 3 -> xor = 0
//                n % 4 == 1 -> xor = 1
//                n % 4 == 2 -> xor = n+1

// Q : for a range [l...r], find xor of all numbers from l to r; 1 <= l <= r <= 1e18
// Sol -> find xor in range [1..l-1], [1..r], then do xor(1..l-1)^xor(1..r) 

// Q : max length subarray whose xor is k. (imp)
// Sol -> thing about simething like prefix xor + hashing

// Power Set: Q - 1097B (practice - codeforces)
// Power Set is a better method to create all the subsequences of a given array or string rather than recursion.
int n; // size of array whose power set is being constructed.
int arr[n]; // array whose power set is being constructed
vector<vector<int>> power_set;
for(int num = 0; num < (1<<n); nums++){
    vector<int> s;
    for(int i = 0; i < n; i++){
        if(num & (1 << i)) s.push_back(arr[i]);
    }
    power_set.push_back(s);
}
// the matrix power_set cantains all the subsequences of given array of size n.

/*.................................................................................................................................................................*/

// combinatorics:

// permutations and combinations:
// calculate nCr for given calue if n and r, and T testcases.
// constraints -> 1 <= T <= 1e5
//             -> 1 <= r <= n <= 1e5
// for brute force the complexity will be (T*N) i.e. 1e10 operations -> TLE

// optimised ideology : when we compute n!, we actually compute factorials for all the number from 1 to n-1 and then find n!.
// so to reduce repetative calculations -> we can find factorial of n = 1e5, and in this way we have factorial of all numbers 
// 1 to n = 1e5, and as we have stored these values calculations will have access to factorial values in O(1) and hence we have 
// optimised the time complexity by doing pre-computation.

// if we have to find nCr for just a single value of n -> use the following function:
#define int long long // defines every int data type to long long.
int mod = 1e9+7;

// modular exponentiation: -> computing (base)^n
int power(int base, int n){ // T.C : O(log(base))
    int ans = 1;
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

void pre_computation(vector<int> &fact){ // T.C : O(n)
    fact[0] = 1;
    for(int i = 1; i<=fact.size(); i++){
        fact[i] = (i * fact[i-1]) % mod;
    }
}

int nCr(int n, int r, vector<int> &fact){ // T.C : O(2 * logn).
    return (fact[n] *(power(fact[r], mod-2) * power(fact[n-r], mod-2)) % mod) % mod;
}

//int main() -> signed main()
signed main(){ // converted int to signed as we have declared int as long long.
    int N = 1e5;
    vector<int> fact(N, 0);
    pre_computation(fact);

    int n, r;
    cin >> n >> r;
    cout << nCr(n, r, fact); 
    return 0;
}

                                            ------------------------------------------------------------------------

// if we have multiple test cases, then pre_computating the power(fact(n)) will be benificail, hence for multiple test cases use following function:
#define int long long // defines every int data type to long long.
int mod = 1e9+7;

// modular exponentiation: -> computing (base)^n
int power(int base, int n){ // T.C : O(log(base))
    int ans = 1;
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

void pre_computation(vector<int> &fact, vector<int> &powerr){ // T.C : O(nlogn)
    fact[0] = 1;
    powerr[0] = 1;
    for(int i = 1; i<=fact.size(); i++){
        fact[i] = (i * fact[i-1]) % mod;
        powerr[i] = power(fact[i], mod-2);
    }
}

int nCr(int n, int r, vector<int> &fact){ // T.C : O(1).
    return (fact[n] *(powerr[r] * powerr[n-r]) % mod) % mod;
}

// int main() -> signed main()
signed main(){ // converted int to signed as we have declared int as long long.
    int N = 1e5;
    vector<int> fact(N, 0);
    vector<int> powerr(N, 0);
    pre_computation(fact, powerr);

    int n, r;
    cin >> n >> r;
    cout << nCr(n, r, fact); 
    return 0;
}

                                            ------------------------------------------------------------------------

// Questions for practice:
// Q : different ways to represent N as a sum of k non-zero numbers.
// example testcase : n = 5 and k = 3 
// -> 3,1,1
// -> 1,3,1
// -> 1,1,3
// -> 2,2,1
// -> 2,1,2
// -> 1,1,2
// Sol : nCr ; here n->(n-1) & r->(k-1).

// Q : grid unique paths - Leetcode : can be done using combinatorics

// Q : k-special cell -> hackerearth (imp)

/*.................................................................................................................................................................*/





/*.................................................................................................................................................................*/
