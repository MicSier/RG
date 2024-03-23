#include <iostream>
#include <symbolicc++.h>
#include <cmath>

using namespace std;

Symbolic slad1(Symbolic x, int n)
{
	Symbolic s = 0;

	for (int i = 0; i < n; i++)
	{
		s = s + x(i);
	}
	return s;
}

Symbolic product(Symbolic a, Symbolic b, int n)
{
	Symbolic p("p", n);

	for (int i = 0; i < n; i++)
	{
		p(i) = a(i) * b(i);
	}
	return p;
}

Symbolic exps(Symbolic x, int n)
{
	Symbolic p("p", n);

	for (int i = 0; i < n; i++)
	{
		p(i) = exp(x(i));
	}
	return p;
}

double log_base(double base, double value)
{
	return log(value) / log(base);
}

double ita(double j1)
{
	return 1.0 / 2.0 * log_base(2.0, 0.25 * exp(3.0 * j1) + 0.75 * exp(-j1)) + 1.0 / 2.0 * log_base(2.0, 0.75 * exp(j1) + 0.25 * exp(-3.0 * j1));
}
double itj(double j1)
{
	return 1.0 / 2.0 * log_base(2.0, 0.25 * exp(3.0 * j1) + 0.75 * exp(-j1)) - 1.0 / 2.0 * log_base(2.0, 0.75 * exp(j1) + 0.25 * exp(-3.0 * j1));
}

double em(double t0)
{
	double nj = 1.0 * t0, nA0 = 0.0, nJ;
	double em = 0.0;
	for (int j = 0; j < 10; j++)
	{
		nA0 = ita(nj);
		nJ = itj(nj);
		em += nA0 / (pow(2.0, j));
		nj = nJ;
	}

	return em;
}

int main()
{
	int L = 4, noc = pow(2, L);
	Symbolic s("s", L, noc), S1("S1"), S2("S2"), P0, e4("e4", 1, noc), h0, j1("j1"), Z, fh, fh0, Z0, Z1, a("a"), b("b"), c("c"), d("d"), Lhs, Rhs, h("h"), t("t");


	for (int i = 0; i < noc; i++)
	{
		e4(0, i) = 1;
	}

	for (int i = 0; i < noc; i++)
	{
		int temp = i;
		for (int j = 0; j < L; j++)
		{
			s(L - j - 1, i) = 2 * ((temp + 1) % 2) - 1;
			temp = temp / 2;
		}
	}

	P0 = product(e4.row(0) + s.row(0) * S1, (e4.row(0) + s.row(3) * S2), noc);
	P0 = P0 / pow(2, L);

	h0 = 0;
	for (int i = 0; i < 3; i++)
	{
		h0 += j1 * product(s.row(i), s.row(i + 1), noc);
	}

	Z = slad1(product(P0, exps(h0, noc), noc), noc);
	Z0 = Z[S1 == 0, S2 == 0];
	Z1 = (Z - Z0)[S1 == 1, S2 == 1];
	//std::cout<<Z-Z0-Z1*S1*S2<<std::endl;
	//std::cout<<Z<<std::endl;
	Lhs = a + b * S1 * S2;
	Rhs = log(2, c + d * S1 * S2);

	Equation eq1 = Lhs[S1 == 1, S2 == 1] == Rhs[S1 == 1, S2 == 1];
	Equation eq2 = Lhs[S1 == -1, S2 == 1] == Rhs[S1 == -1, S2 == 1];

	Equations eqs = (eq1, eq2);
	PatternMatches res = solve(eqs, (a, b));
	Symbolic iter[2];
	int i = 0;
	for (auto x : res)
	{
		for (auto x1 : x)
		{
			iter[i] = x1.rhs[c == Z0, d == Z1];
			i++;
		}

	}
	std::cout << iter[0] << std::endl;
	std::cout << iter[1] << std::endl;
	
	double dh = 0.00001, t0 = 0.2, ds = 0.05;
	for (int i = 0; i < 20; i++)
	{
		double emp = em(t0 - dh);
		double emc = em(t0);
		double emn = em(t0 + dh);
		double ew = (emn-emp)/(2*dh);
		double cw = (emn - 2*emc + emp) / (dh * dh);
		cout << 1.0/t0 << " " << emc << " " << ew << " " << cw << endl;
		t0 += ds;
	}


	return 0;
}
