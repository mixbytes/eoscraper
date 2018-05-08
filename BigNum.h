#pragma once
#include <algorithm>
#include <array>
#include <type_traits>
#include <cstddef>
#include <assert.h>


namespace BigNum
{
	namespace Core
	{
		template<typename T, int N>
		void zero(T (&to)[N]) noexcept
		{
			for (int i = 0; i < N; i++)
				to[i] = 0;
		}

		template<typename T, int N, int M>
		void copy(T (&to)[N], const T (&from)[M]) noexcept
		{
			for (int i = 0; i < N; i++)
				to[i] = (i < M) ? from[i] : 0;
		}

		template<typename T>
		int cmp(const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;
			const std::ptrdiff_t vmax = (v1len > v2len)? v1len : v2len;

			for(std::ptrdiff_t i = vmax; i-->0;)
			{
				const T v1 = (i < v1len) ? v1begin[i] : 0;
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				if (v1 != v2)
					return v1 > v2 ? 1 : -1;
			}
			return 0;
		}

        template<typename T, int N, int M>
      	int cmp(const T(&v1)[N], const T(&v2)[M]) noexcept
		{
			const int v1len = N;
			const int v2len = M;
			const int vmax = (v1len > v2len)? v1len : v2len;

			for(int i = vmax; i-->0;)
			{
				const T v1v = (i < v1len) ? v1[i] : 0;
				const T v2v = (i < v2len) ? v2[i] : 0;
				if (v1v != v2v)
					return v1v > v2v ? 1 : -1;
			}
			return 0;
		}

		template<typename T, int N>
		int cmp(const T(&v1)[N], const T v2) noexcept
		{
			const int v1len = N;
			const int v2len = 1;
			const int vmax = v1len;

			for (int i = vmax; i-->0;)
			{
				const T v1v = (i < v1len) ? v1[i] : 0;
				const T v2v = (i < v2len) ? v2 : 0;
				if (v1v != v2v)
					return v1v > v2v ? 1 : -1;
			}
			return 0;
		}

		template<typename T>
		T adc(T& ret, const T a, const T b) noexcept
		{
			ret = a + b;
			return (ret < a) ? 1 : 0;
		}

		unsigned long long int adc(unsigned long long int& ret, const unsigned long long int v1, const unsigned long long int v2, const unsigned long long int carry) noexcept
		{
			unsigned long long int retCarry = adc(ret, v1, v2);
			retCarry = adc(ret, ret, carry) + retCarry;
			return retCarry;
		}

		template <typename T>
		T adc(T& ret, const T v1, const T v2, const T carry, typename std::enable_if<(sizeof(unsigned long long int) > sizeof(T))>::type* = 0) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			const unsigned long long int kWordMaskBits = (1ull << (kWordSizeBits)) - 1ull;
			static_assert(sizeof(unsigned long long int) > sizeof(T), "cant detect overflow");
			
            const unsigned long long int sum = (unsigned long long int)(v1) + v2 + carry;
            ret = sum & kWordMaskBits;
            return T(sum >> kWordSizeBits);
		}

		template <typename T>
		T add(T* rezbegin, T* rezend, const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v1 = (i < v1len) ? v1begin[i] : 0;
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}
		
		template<typename T, int N, int M, int K>
		T add(T (&rez)[N], const T (&v1)[M], const T (&v2)[K]) noexcept
		{
			T carry = 0;
			for (int i = 0; i < N; i++)
			{
				const T v1v = (i < M) ? v1[i] : 0;
				const T v2v = (i < K) ? v2[i] : 0;
				carry = adc(rez[i], v1v, v2v, carry);
			}
			return carry;
		}

		template<typename T>
		T add(T* rezbegin, T* rezend, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = adc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		template<typename T>
		T sbc(T& ret, const T a, const T b) noexcept
		{
			ret = a - b;
			return (ret > a) ? 1 : 0;
		}

		unsigned long long int sbc(unsigned long long int& ret, const unsigned long long int v1, const unsigned long long int v2, const unsigned long long int carry) noexcept
		{
			unsigned long long int retCarry = sbc(ret, v1, v2);
			retCarry = sbc(ret, ret, carry) + retCarry;
			return retCarry;
		}

		template<typename T>
		T sbc(T& ret, const T v1, const T v2, const T carry, typename std::enable_if<(sizeof(unsigned long long int) > sizeof(T))>::type* = 0) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			const unsigned long long int kWordMaskBits = (1ull << (kWordSizeBits)) - 1ull;
			static_assert(sizeof(unsigned long long int) > sizeof(T), "cant detect overflow");
			
			const unsigned long long int sum = (unsigned long long int)(v1) - v2 - carry;
			ret = sum & kWordMaskBits;
			return (sum >> kWordSizeBits) ? 1u : 0u;
		}

		template<typename T>
		T sub(T* rezbegin, T* rezend, const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v1 = (i < v1len) ? v1begin[i] : 0;
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], v1, v2, carry);
			}
			return carry;
		}

		template<typename T, int N, int M, int K>
		T sub(T (&rez)[N], const T (&v1)[M], const T (&v2)[K]) noexcept
		{
			T carry = 0;
			for (int i = 0; i < N; i++)
			{
				const T v1v = (i < M) ? v1[i] : 0;
				const T v2v = (i < K) ? v2[i] : 0;
				carry = sbc(rez[i], v1v, v2v, carry);
			}
			return carry;
		}

		template<typename T>
		T sub(T* rezbegin, T* rezend, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t rezlen = rezend - rezbegin;
			const std::ptrdiff_t v2len = v2end - v2begin;

			T carry = 0;
			for (std::ptrdiff_t i = 0; i < rezlen; i++)
			{
				const T v2 = (i < v2len) ? v2begin[i] : 0;
				carry = sbc(rezbegin[i], rezbegin[i], v2, carry);
			}
			return carry;
		}

		unsigned long long int _mulhi(unsigned long long int& ret, const unsigned long long int v1, const unsigned long long int v2)
		{
			const unsigned int kWordSizeBits = sizeof(unsigned long long int) * 8;
			const unsigned int kHalfWordSizeBits = kWordSizeBits / 2;
			const unsigned long long int kHalfWordMaskBits = (1ull << (kHalfWordSizeBits)) - 1ull;
			
			unsigned long long int v1l = v1 & kHalfWordMaskBits;
			unsigned long long int v1h = v1 >> kHalfWordSizeBits;
			unsigned long long int v2l = v2 & kHalfWordMaskBits;
			unsigned long long int v2h = v2 >> kHalfWordSizeBits;

			unsigned long long int z0 = v1l * v2l;
			unsigned long long int z1 = v1l * v2h;
			unsigned long long int z2 = v1h * v2l;
			unsigned long long int z3 = v1h * v2h;

			unsigned long long int carry1 = adc(ret, z0, z1 << kHalfWordSizeBits);
			unsigned long long int carry2 = adc(ret, ret, z2 << kHalfWordSizeBits);
			unsigned long long int hi = 0;
			adc(hi, z1 >> kHalfWordSizeBits, z2 >> kHalfWordSizeBits, carry1);
			adc(hi, hi, z3, carry2);
			return hi;
		}

		template<typename T>
		inline T madd(T& ret, const T v1, const T v2, const T carry, typename std::enable_if<sizeof(unsigned long long int) == sizeof(T)>::type* = 0)
		{
			T ml = 0;
			T mh = _mulhi(ml, v1, v2);
			T mc = adc(ret, ret, ml);
			mc += adc(ret, ret, carry);
			return mh + mc;	
		}

		template<typename T>
		inline T madd(T& ret, const T v1, const T v2, const T carry, typename std::enable_if<(sizeof(unsigned long long int) > sizeof(T))>::type* = 0)
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;

			unsigned long long int m = (unsigned long long int)(ret) + (unsigned long long int)(v1) * (unsigned long long int)(v2) + carry;
			ret = T(m);
			return T(m >> kWordSizeBits);
		}

		template<typename T>
		void mul(T* rezbegin, T* rezend, const T* v1begin, const T* v1end, const T* v2begin, const T* v2end) noexcept
		{
			const std::ptrdiff_t v1len = v1end - v1begin;
			const std::ptrdiff_t v2len = v2end - v2begin;
			const std::ptrdiff_t rezlen = rezend - rezbegin;

			assert((rezlen) >= (v1len) + (v2len));
			assert(std::any_of(rezbegin, rezend, [&](T w) { return w > 0u; }) == false);

			const ptrdiff_t v1last = std::min(rezlen, v1len);
			for (std::ptrdiff_t v1it = 0; v1it < v1last; ++v1it)
			{
				const ptrdiff_t v2last = std::min(v2len, rezlen - v1it);
				T carry = 0;
				for (std::ptrdiff_t v2it = 0; v2it < v2last; ++v2it)
					carry = madd(rezbegin[v1it + v2it], v1begin[v1it], v2begin[v2it], carry);
				rezbegin[v1it + v2last] = carry;
			}
		}
		
		template<typename T>
		void mul(T(&rez)[1], const T(&v1)[1], const T(&v2)[1]) noexcept
		{
			rez[0] = v1[0] * v2[0];
		}

		template<typename T, int N, int M, int K>
		void mul(T (&rez)[N], const T (&v1)[M], const T (&v2)[K], typename std::enable_if<N >= M + K>::type* = 0) noexcept
		{
			assert(std::any_of(std::begin(rez), std::end(rez), [&](T w) { return w > 0; }) == false);

			for (int v1it = 0; v1it < M; ++v1it)
			{
				T carry = 0;
				for (int v2it = 0; v2it < K; ++v2it)
					carry = madd(rez[v1it + v2it], v1[v1it], v2[v2it], carry);
				rez[v1it + K] = carry;
			}
		}

		template<typename T, int N, int M, int K>
		void mul(T(&rez)[N], const T(&v1)[M], const T(&v2)[K], typename std::enable_if<N < M + K>::type* = 0) noexcept
		{
			assert(std::any_of(std::begin(rez), std::end(rez), [&](T w) { return w > 0; }) == false);

			const int v1last = std::min(N, M);
			for (int v1it = 0; v1it < v1last; ++v1it)
			{
				const int v2last = std::min(K, N - v1it);
				T carry = 0;
				for (int v2it = 0; v2it < v2last; ++v2it)
					carry = madd(rez[v1it + v2it], v1[v1it], v2[v2it], carry);
				rez[v1it + v2last] = carry;
			}
		}
		
		template<typename T, unsigned int N>
		void shl(T (&rez)[N], const T (&v)[N], unsigned int count) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			
			assert(count < N*sizeof(T)*8);
			const unsigned int bytes = count / kWordSizeBits;
			const unsigned int bits = count % kWordSizeBits;
			if (bits == 0)
			{
				for (unsigned int i = N; i --> bytes;)
					rez[i] = v[i - bytes];
			}
			else
			{
				unsigned int i = 0;
				for (; i < bytes; i++)
					rez[i] = 0;
				rez[i] = (v[i - bytes] << (bits));
				i++;
				for (; i < N; i++)
					rez[i] = (v[i - bytes] << (bits)) | (v[i - bytes - 1] >> (sizeof(T) * 8 - bits));
			}
		}
		
		template<typename T, unsigned int N>
		void shr(T (&rez)[N], const T (&v)[N], unsigned int count) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			assert(count < N*sizeof(T)*8);
			const unsigned int bytes = count / kWordSizeBits;
			const unsigned int bits = count % kWordSizeBits;
			if (bits == 0)
			{
				for (unsigned int i = 0; i < N - bytes; i++)
					rez[i] = v[i + bytes];
			}
			else
			{
				unsigned int i = 0;
				for (; i < N - bytes - 1; i++)
					rez[i] = (v[i + bytes + 1] << (sizeof(T) * 8 - bits)) | (v[i + bytes] >> bits);
				rez[i] = (v[i + bytes] >> bits);
				i++;
				for (; i < N; i++)
					rez[i] = 0;
			}
		}

		template<typename T, int N>
		void calcShiftedD(T(&ret)[sizeof(T) * 8][N + 1], const T(&d)[N]) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			
			for(unsigned i = 0; i < N; i++)
				ret[0][i] = d[i];
			ret[0][N] = 0;

			for (unsigned i = 1; i < kWordSizeBits; i++)
				shl(ret[i], ret[0], i);
		}

		template<typename T, int N, int M>
		void div(T(&q)[N], T(&n)[N], const T(&d)[M]) noexcept
		{
			if (cmp(n,d) < 0)
			{
				zero(q);
				return;
			}
			
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			
			T shiftedD[kWordSizeBits][M + 1];
			calcShiftedD(shiftedD, d);

			T* rbeg = std::end(n);
			T* rend = std::end(n);

			for (int i = N; i-->0;)
			{
				rbeg--;
				T val = 0;
				if (Core::cmp(rbeg, rend, std::begin(d), std::end(d)) >= 0)
				{
					for (int j = kWordSizeBits; j-- > 0;)
					{
						if (Core::cmp(rbeg, rend, std::begin(shiftedD[j]), std::end(shiftedD[j])) >= 0)
						{
							Core::sub(rbeg, rend, std::begin(shiftedD[j]), std::end(shiftedD[j]));
							val += T(1) << j;
						}
					}
				}
				q[i] = val;
			}
		}

		template<typename T, int N>
		void div(T(&q)[N], T(&n)[N], const T d, typename std::enable_if<(sizeof(unsigned long long int) > sizeof(T))>::type* = 0) noexcept
		{
			static_assert(sizeof(unsigned long long) > sizeof(T), "T have to be smaller size than unsigned long long");
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			
			unsigned long long rest = 0;
			for (int i = N; i--> 0;)
			{
				rest = rest << kWordSizeBits;
				rest += n[i];
				n[i] = 0;
				q[i] = T(rest / d);
				rest = rest % d;
			}
			n[0] = (T)rest;
		}
		
		template<typename T, int N>
		void div(T(&q)[N], T(&n)[N], const T d, typename std::enable_if<(sizeof(unsigned long long int) == sizeof(T))>::type* = 0) noexcept
		{
			const T converted_d[] = { d };
			div(q, n, converted_d);
		}

		template<typename T, int N, int M>
		void mod(T(&n)[N], const T(&d)[M]) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			
			T shiftedD[kWordSizeBits][M + 1];
			calcShiftedD(shiftedD, d);

			T* rbeg = std::end(n);
			T* rend = std::end(n);

			for (int i = N; i-->0;)
			{
				rbeg--;
				if (Core::cmp(rbeg, rend, std::begin(d), std::end(d)) >= 0)
				{
					for (int j = kWordSizeBits; j-- > 0;)
					{
						if (Core::cmp(rbeg, rend, std::begin(shiftedD[j]), std::end(shiftedD[j])) >= 0)
							Core::sub(rbeg, rend, std::begin(shiftedD[j]), std::end(shiftedD[j]));
					}
				}
			}
		}

		template<typename T, int N>
		T mod(const T(&n)[N], const T d, typename std::enable_if<(sizeof(unsigned long long int) > sizeof(T))>::type* = 0) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			static_assert(sizeof(unsigned long long) > sizeof(T), "T have to be smaller size than unsigned long long");
			unsigned long long rest = 0;
			for (int i = N; i--> 0;)
			{
				rest = rest << kWordSizeBits;
				rest += n[i];
				rest = rest % d;
			}
			return (T)rest;
		}
		
		template<typename T, int N>
		T mod(const T(&n)[N], const T d, typename std::enable_if<(sizeof(unsigned long long int) == sizeof(T))>::type* = 0) noexcept
		{
			const T tmp_d[] = { d };
			T tmp_n[N];
			copy(tmp_n, n);
			mod(tmp_n, tmp_d);
			return tmp_n[0];
		}

		template<typename T, int N>
		void gcd(T(&rez)[N], const T(&v1)[N], const T(&v2)[N]) noexcept
		{
			int cmprez = cmp(v1, v2);
			if (cmprez > 0)
				copy(rez, v1);
			else
				copy(rez, v2);
			
			if (cmprez == 0)
				return;

			T b[N];
			if (cmprez > 0)
				copy(b, v2);
			else
				copy(b, v1);

			while (cmp(b, T(0)))
			{
				T t[N], m[N];
				copy(t, b);
				copy(m, rez);
				mod(m, b);
				copy(b, m);
				copy(rez, t);
			}
		}
		
		template<typename T, int N>
		void lcm(T(&rez)[2*N], const T(&v1)[N], const T(&v2)[N]) noexcept
		{
			T m[2 * N];
			zero(m);
			mul(m, v1, v2);
			T g[N];
			gcd(g, v1, v2);
			div(rez, m, g);
		}
		
		template<typename T>
		const T modInv(const T& a, const T& n) noexcept
		{
			T t(0), newt(1);
			T r = n, newr = a;
			while (newr != T(0))
			{
				T quotient = r / newr;
				{
					//(t, newt) := (newt, t - quotient * newt) 
					T tmp = newt;
					newt = t - quotient * newt;
					t = tmp;
				}
				{
					//(r, newr) := (newr, r - quotient * newr)
					T tmp = newr;
					newr = r - quotient * newr;
					r = tmp;
				}
			}
			if (r > T(1)) return T(0);
			if (t < 0 || t > n) t = t + n;
			return t;
		}
		
		template<typename T, int N>
		void modInv(T(&rez)[N], const T(&a)[N], const T(&n)[N]) noexcept
		{
		    T t[N], newt[N];
		    zero(t);
		    zero(newt);
		    newt[0] = T(1);
		    
		    T r[N], newr[N];
		    copy(r, n);
		    copy(newr, a);
		    
		    while (cmp(newr, T(0)) != 0)
		    {
		        T quotient[N];
		        T tmpr[N];
		        copy(tmpr, r);
		        div(quotient, tmpr, newr);
		        {
		        	//(t, newt) := (newt, t - quotient * newt) 
		        	T tmp[N];
		        	copy(tmp, newt);
		        	T qn[2*N];
		        	zero(qn);
		        	mul(qn, quotient, newt);
		        	sub(newt, t, qn);
		        	copy(t,tmp);
		        }
		        {
		        	//(r, newr) := (newr, r - quotient * newr)
		        	T tmp[N];
		        	copy(tmp, newr);
		        	T qn[2*N];
		        	zero(qn);
		        	mul(qn, quotient, newr);
		        	sub(newr, r, qn);
		        	copy(r,tmp);
		        }
		    }
		    if (cmp(r, T(1)) > 0)
		    {
		    	zero(rez);
		    	return;
		    }
			if (cmp(t, n) > 0)
			{
				add(t, t, n);
			}
		    copy(rez, t);
		}
		
		template<typename T, int N, int M, unsigned int K>
		void modExp(T (&result)[N], const T (&base)[N], const T (&exp)[K], const T (&modulo)[M]) noexcept
		{
			const unsigned int kWordSizeBits = sizeof(T) * 8;
			
			zero(result);
			if (cmp(modulo,T(1)) == 0)
				return;
	
			result[0] = T(1);

			unsigned int realExpSize = K;
			// skip leading zeros
			while (realExpSize > 0 && exp[realExpSize - 1] == 0)
				realExpSize--;
	
			for (unsigned int i = realExpSize; i-->0;)
			{
				for (unsigned int j = (kWordSizeBits); j-->0;)
				{
					T tmp[2*N], q[2*N];
					zero(tmp); zero(q);
					mul(tmp, result, result);
					div(q, tmp, modulo);
					copy(result, tmp);

					if ((exp[i] >> j) & T(1))
					{
						zero(tmp); zero(q);
						mul(tmp, result, base);
						div(q, tmp, modulo);
						copy(result, tmp);
					}
				}
			}
		}
	}
	
	typedef unsigned long long int Word;
	const unsigned int kWordSizeBits = sizeof(Word) * 8;
	
	template<int N>
	class Num
	{
	public:
        const uint64_t kWordSizeBits ;
        const uint64_t Size ;

        Word data[(N + (sizeof(Word) * 8)-1) / (sizeof(Word) * 8)];

        Num() : Size((N + (kWordSizeBits)-1) / (kWordSizeBits)), kWordSizeBits(sizeof(Word) * 8)
		{
			for (unsigned int i = 0; i < Size; i++)
				data[i] = 0;
		}

		template<typename T>
		Num(T num, typename std::enable_if<(sizeof(T) > sizeof(Word))>::type* = 0) : Size((N + (kWordSizeBits)-1) / (kWordSizeBits)), kWordSizeBits(sizeof(Word) * 8)
		{
			for (unsigned int i = 0; i < Size; i++) {
				data[i] = Word(num);
				num >>= kWordSizeBits;
			}
		}
		
		template<typename T>
		Num(T num, typename std::enable_if<(sizeof(T) <= sizeof(Word))>::type* = 0) : Size((N + (kWordSizeBits)-1) / (kWordSizeBits)), kWordSizeBits(sizeof(Word) * 8)
		{
			data[0] = Word(num);
			for (unsigned int i = 1; i < Size; i++)
				data[i] = 0;
		}

		template<int M>
		Num(const Num<M>& other)  : Size((N + (kWordSizeBits)-1) / (kWordSizeBits)), kWordSizeBits(sizeof(Word) * 8)
		{
			*this = other;
		}

		template<int M>
		Num(Num<M>&& other)  : Size((N + (kWordSizeBits)-1) / (kWordSizeBits)), kWordSizeBits(sizeof(Word) * 8)
		{
			*this = other;
		}

		Num<N>& operator = (const Num<N>& other) noexcept
		{
			if (this != &other)
			{
				for (unsigned int i = 0; i < Size; i++)
					data[i] = other[i];
			}
			return *this;
		}

		template<int M>
		Num<N>& operator = (const Num<M>& other) noexcept
		{
			for (unsigned int i = 0; i < std::min(Size, other.Size); i++)
				data[i] = other[i];

			for (unsigned int i = std::min(Size, other.Size); i < Size; i++)
				data[i] = 0;
				
			//assert(std::none_of(&other[std::min<int>(Size, Num<M>::Size)],&other[Num<M>::Size], [](unsigned char x){ return x != 0; }));

			return *this;
		}

		Word& operator[](unsigned int i) noexcept { return data[i]; }
		const Word& operator[](unsigned int i) const noexcept { return data[i]; }

		bool bit(unsigned int i) const noexcept
		{
			assert(i >= 0 && i < N);
			bool ret = false;
			if (i >= 0 && i < N)
			{
				const unsigned int byte = i / kWordSizeBits;
				const unsigned int bit = i % kWordSizeBits;
				ret = ((*this)[byte] & (Word(1) << bit)) != 0;
			}
			return ret;
		}

		Word* begin() noexcept { return &data[0]; }
		Word* end() noexcept { return &data[Size]; }

		const Word* begin() const noexcept { return &data[0]; }
		const Word* end() const noexcept { return &data[Size]; }
	};

	template<int N, int M>
	const Num<N> operator + (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> ret;
		Core::add(ret.data, v1.data, v2.data);
		return ret;
	}

	template<int N>
	const Num<N> operator + (const Num<N>& v1, Word v2) noexcept
	{
		Num<N> ret;
		const Word a2[] = {v2};
		Core::add(ret.data, v1.data, a2);
		return ret;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v) noexcept
	{
		return Num<N>(0) - v;
	}

	template<int N, int M>
	const Num<N> operator - (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> ret;
		Core::sub(ret.data, v1.data, v2.data);
		return ret;
	}

	template<int N>
	const Num<N> operator - (const Num<N>& v1, Word v2) noexcept
	{
		Num<N> ret;
		const Word a2[] = { v2 };
		Core::sub(ret.data, v1.data, a2);
		return ret;
	}

	template<int N, int M>
	const Num<N + M> operator * (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N + M> rez(0);
		Core::mul(rez.data, v1.data, v2.data);
		return rez;
	}

	template<int N>
	const Num<N + kWordSizeBits> operator * (const Num<N>& v1, Word v2) noexcept
	{
		Num<N + kWordSizeBits> rez(0);
		const Word a2[] = { v2 };
		Core::mul(rez.data, v1.data, a2);
		return rez;
	}

	template<int N, int M>
	const bool operator == (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) == 0;
	}

	template<int N>
	const bool operator == (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) == 0;
	}

	template<int N, int M>
	const bool operator != (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return !(v1 == v2);
	}
	
	template<int N>
	const bool operator != (const Num<N>& v1, Word v2) noexcept
	{
		return !(v1 == v2);
	}

	template<int N, int M>
	const bool operator > (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) > 0;
	}

	template<int N>
	const bool operator > (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) > 0;
	}

	template<int N, int M>
	const bool operator >= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) >= 0;
	}

	template<int N>
	const bool operator >= (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) >= 0;
	}

	template<int N, int M>
	const bool operator < (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) < 0;
	}

	template<int N>
	const bool operator < (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) < 0;
	}

	template<int N, int M>
	const bool operator <= (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return Core::cmp(v1.data, v2.data) <= 0;
	}

	template<int N>
	const bool operator <= (const Num<N>& v1, Word v2) noexcept
	{
		const Word a2[] = { v2 };
		return Core::cmp(v1.data, a2) <= 0;
	}

	template<int N>
	std::ostream & operator << (std::ostream & out, const Num<N>& num) {
		std::string ret;
		BigNum::Num<N> Q = num;

		while (Q > 0)
		{
			auto R = Q % 10;
			assert(R < 10);
			Q = Q / 10;
			ret.push_back(char(R + '0'));
		}
		std::reverse(ret.begin(), ret.end());
		if (ret.empty())
			ret = "0";
		return (out << ret);
	}

	template<int N>
	const Num<N> operator << (const Num<N>& v1, unsigned int bits) noexcept
	{
		Num<N> ret;
		if (bits < N)
			Core::shl(ret.data, v1.data, bits);
		return ret;
	}

	template<int N>
	const Num<N> operator >> (const Num<N>& v1, unsigned int bits) noexcept
	{
		Num<N> ret;
		if (bits < N)
			Core::shr(ret.data, v1.data, bits);
		return ret;
	}

	template<int N, int M>
	const Num<N> operator / (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		if (v1 < v2)
			return Num<N>(0);

		Num<N> q, r = v1;
		Core::div(q.data, r.data, v2.data);
		return q;
	}

	template<int N>
	const Num<N> operator / (const Num<N>& v1, const Word v2) noexcept
	{
		Num<N> q, r = v1;
		Core::div(q.data, r.data, v2);
		return q;
	}

	template<int N, int M>
	const Num<M> operator % (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		if (v1 < v2)
			return Num<M>(v1);

		Num<N> r = v1;
		Core::mod(r.data, v2.data);
		return r;
	}

	template<int N>
	const Word operator % (const Num<N>& v1, const Word v2) noexcept
	{
		return Core::mod(v1.data, v2);
	}

	template<int N, int M>
	const Num<N> operator & (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		Num<N> rez(0);

		for (int i = 0; i < v1.Size; i++)
		{
			rez[i] = v1[i] & ((i < M) ? v2[i] : 0);
		}

		return rez;
	}

	template<int N, int M>
	const Num<N> operator & (const Num<N>& v1, Num<M>&& v2) noexcept
	{
		Num<N> rez(0);

		for (int i = 0; i < v1.Size; i++)
		{
			rez[i] = v1[i] & ((i < M) ? v2[i] : 0);
		}

		return rez;
	}

	template<int N>
	const Num<N> operator & (const Num<N>& v1, const Word v2) noexcept
	{
		Num<N> rez(v1[0] & v2);
		return rez;
	}

	template<int N>
	const Num<N> operator | (const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<N> rez;

		for (unsigned int i = 0; i < v1.Size; i++)
		{
			rez[i] = v1[i] | v2[i];
		}

		return rez;
	}

	template<int N>
	const Num<N> operator | (const Num<N>& v1, const Word v2) noexcept
	{
		Num<N> rez = v1;
		rez[0] = v1[0] | v2;
		return rez;
	}

	template<int N>
	const Num<N> operator ~ (const Num<N>& v) noexcept
	{
		Num<N> rez;

		for (unsigned int i = 0; i < v.Size; i++)
		{
			rez[i] = ~v[i];
		}

		return rez;
	}

	template<int N, int M>
	bool operator && (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return v1 != 0 && v2 != 0;
	}

	template<int N, int M>
	bool operator || (const Num<N>& v1, const Num<M>& v2) noexcept
	{
		return v1 != 0 || v2 != 0;
	}

	template<int N>
	bool operator !(const Num<N>& v) noexcept
	{
		return v == 0;
	}
	
	template<int N>
	const Num<N> gcd(const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<N> ret;
		Core::gcd(ret.data, v1.data, v2.data);
		return ret;
	}

	template<int N>
	const Num<2*N> lcm(const Num<N>& v1, const Num<N>& v2) noexcept
	{
		Num<2 * N> ret;
		Core::lcm(ret.data, v1.data, v2.data);
		return ret;
	}

	template<int N>
	const Num<N> modInv(const Num<N>& a, const Num<N>& n) noexcept
	{
	    Num<N> rez;
	    Core::modInv(rez.data, a.data, n.data);
	    return rez;
	}

	template<typename T>
	const T modInv(const T& a, const T& n) noexcept
	{
		return Core::modInv(a, n);
	}

	template<int N>
	const Num<N> modExp(const Num<N>& base, const Num<N>& exp, const Num<N>& modulo) noexcept
	{
		Num<N> result;
		Core::modExp(result.data, base.data, exp.data, modulo.data);
		return result;
	}

};