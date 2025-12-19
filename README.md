# Numerical Methods
### Introduction


### Project Structure

```
README.md/
├── Project Introduction
├── README.md Structure
└── Table of Contents
├── Solution of Non-Linear Equations
    ├── Bisection
    ├── False Position
    ├── Secant
    ├── Newton Raphson

├── Solution of Linear Equations
    ├── Gauss Elimination
    ├── Gauss Jordan Elimination
    ├── LU Decomposition
    ├── Matrix Inversion

├── Solution of Differential Equations
    ├── Runge Kutta 4th Order

├── Numerical Integration
    ├── Simpson's 1/3 Rule
    └── Simpson's 3/8 Rule

├── Interpolation Methods
    ├── Newton Forward Interpolation
    ├── Newton Backward Interpolation
    ├── Newton Divided Difference Interpolation

├── Numerical Differentiation
    ├── Differentiation by Forward Interpolation
    ├── Differentiation by Backward Interpolation

├── Curve Fitting & Regression
    ├── Linear Regression
    ├── Polynomial Regression
    ├── Transcendental Regression


```


#### Table of Contents

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Differential Equations](#solution-of-differential-equations)
  - [Runge-Kutta 4th Order Method](#runge-kutta-4th-order-method)
    - [Theory](#runge-kutta-4th-order-theory)
    - [Code](#runge-kutta-4th-order-code)
    - [Input](#runge-kutta-4th-order-input)
    - [Output](#runge-kutta-4th-order-output)

- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule](#simpsons-13-rule)
    - [Theory](#simpsons-13-theory)
    - [Code](#simpsons-13-code)
    - [Input](#simpsons-13-input)
    - [Output](#simpsons-13-output)
  - [Simpson's 3/8 Rule](#simpsons-38-rule)
    - [Theory](#simpsons-38-theory)
    - [Code](#simpsons-38-code)
    - [Input](#simpsons-38-input)
    - [Output](#simpsons-38-output)

- [Interpolation Methods](#interpolation-methods)
  - [Newton Forward Interpolation](#newton-forward-interpolation)
    - [Theory](#newton-forward-theory)
    - [Code](#newton-forward-code)
    - [Input](#newton-forward-input)
    - [Output](#newton-forward-output)
  - [Newton Backward Interpolation](#newton-backward-interpolation)
    - [Theory](#newton-backward-theory)
    - [Code](#newton-backward-code)
    - [Input](#newton-backward-input)
    - [Output](#newton-backward-output)
  - [Newton Divided Difference Interpolation](#newton-divided-difference-interpolation)
    - [Theory](#newton-divided-difference-theory)
    - [Code](#newton-divided-difference-code)
    - [Input](#newton-divided-difference-input)
    - [Output](#newton-divided-difference-output)

- [Numerical Differentiation](#numerical-differentiation)
  - [Differentiation by Forward Interpolation](#differentiation-forward-interpolation)
    - [Theory](#differentiation-forward-theory)
    - [Code](#differentiation-forward-code)
    - [Input](#differentiation-forward-input)
    - [Output](#differentiation-forward-output)
  - [Differentiation by Backward Interpolation](#differentiation-backward-interpolation)
    - [Theory](#differentiation-backward-theory)
    - [Code](#differentiation-backward-code)
    - [Input](#differentiation-backward-input)
    - [Output](#differentiation-backward-output)

- [Curve Fitting / Regression](#curve-fitting--regression)
  - [Linear Regression](#linear-regression)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)
  - [Polynomial Regression](#polynomial-regression)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)
  - [Transcendental Regression](#transcendental-regression)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)


    
   
    
---


### Solution of Non-Linear Equations

### Bisection Method

#### Bisection Theory
Bisection Method

Theory

The Bisection Method is a numerical method which is used to find the real roots of a polynomial function. It works by dividing an interval into two halves and selecting the sub-interval in which the root lies. The method is based on the Intermediate Value Theorem, which guarantees the existence of a root within an interval where the function values at the endpoints have opposite signs. This method gradually reduces the interval containing the root by taking the average (midpoint) of the interval endpoints. Although simple and reliable, the bisection method converges relatively slowly compared to other numerical methods.

Mathematical Principle
Let f(x) be a continuous function defined on a closed interval [a, b], such that:
f(a) * f(b) < 0
Then, according to the Intermediate Value Theorem, there exists at least one value x in (a, b) for which: f(x) = 0

Bisection Method Algorithm
To find the root of a continuous function f(x), follow these steps:

1. Choose two initial points a and b such that a < b and f(a) * f(b) < 0
2. Compute the midpoint:
   t = (a + b) / 2
3. If f(t) = 0, then t is the root
4. Otherwise:
   - If f(t) * f(a) < 0, the root lies in the interval [a, t]
   - Else if f(t) * f(b) < 0, the root lies in the interval [t, b]
5. Repeat the above steps until the interval becomes sufficiently small or the desired accuracy is achieved

The bisection method is an approximation technique that continues dividing the interval until the root is obtained within a set error margin.

Features
- Supports multiple test cases
- Dynamic evaluation of polynomial functions
- Automatic generation of function strings
- Eliminates duplicate roots
- Well-formatted output with precision control

#### Bisection Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> coef;
int deg;

double fx(double x)
{
    double s = 0;
    for (int i = 0; i <= deg; i++)
        s += coef[i] * pow(x, deg - i);
    return s;
}

string funcStr()
{
    stringstream ss;
    ss << "f(x)=";
    bool first = true;

    for (int i = 0; i <= deg; i++)
    {
        double c = coef[i];
        if (fabs(c) < 1e-10) continue;

        if (!first)
            ss << (c > 0 ? "+" : "-");
        else if (c < 0)
            ss << "-";

        double v = fabs(c);
        int p = deg - i;

        if (p == 0)
            ss << v;
        else if (p == 1)
            ss << (fabs(v - 1.0) < 1e-10 ? "x" : to_string(v) + "x");
        else
            ss << (fabs(v - 1.0) < 1e-10 ? "x^" : to_string(v) + "x^") << p;

        first = false;
    }
    return ss.str();
}

double limitCalc()
{
    if (deg < 2) return 10.0;
    double a = coef[0], b = coef[1], c = coef[2];
    if (fabs(a) < 1e-10) return 10.0;
    return sqrt((b / a) * (b / a) - 2 * (c / a));
}

void uniq(vector<double> &r, double e)
{
    sort(r.begin(), r.end());
    r.erase(unique(r.begin(), r.end(),
        [&](double x, double y)
        {
            return fabs(x - y) < e;
        }), r.end());
}

void solve(int tc, ofstream &out)
{
    vector<double> roots;
    double eps = 0.001;
    double lim = limitCalc();

    for (double i = -lim; i <= lim; i += 0.5)
    {
        double a = i, b = i + 0.5;

        if (fx(a) * fx(b) < 0)
        {
            while (fabs(b - a) > eps)
            {
                double m = (a + b) / 2;
                if (fx(a) * fx(m) < 0)
                    b = m;
                else
                    a = m;
            }
            roots.push_back((a + b) / 2);
        }
    }

    uniq(roots, eps);

    out << "TestCase " << tc << "\n";
    out << "Function:" << funcStr() << "\n";
    out << "Degree:" << deg << "\n";
    out << "Tolerance:" << eps << "\n";
    out << "Roots:\n";

    if (roots.empty())
        out << "None\n";
    else
        for (double x : roots)
            out << fixed << setprecision(6) << x << "\n";

    out << "\n";
}

int main()
{
    ifstream in("Input.txt");
    ofstream out("Output.txt");

    int T;
    in >> T;
    out << "BisectionMethod\n";
    out << "Cases:" << T << "\n\n";

    for (int t = 1; t <= T; t++)
    {
        in >> deg;
        coef.assign(deg + 1, 0);

        for (int i = 0; i <= deg; i++)
            in >> coef[i];

        solve(t, out);
    }

    return 0;
}

```

#### Bisection Input
```
3
3
1 -6 11 -6
2
1 -4 4
5
1 -3 -7 27 -10 -8

```

#### Bisection Output
```
BisectionMethod
Cases:3

TestCase:1
Function:f(x)=x^3-6x^2+11x-6
Degree:3
Tolerance:0.001
Roots:
1.000000
2.000000
3.000000

TestCase:2
Function:f(x)=x^2-4x+4
Degree:2
Tolerance:0.001
Roots:
2.000000

TestCase:3
Function:f(x)=x^5-3x^4-7x^3+27x^2-10x-8
Degree:5
Tolerance:0.001
Roots:
-2.000244
0.999756
2.000244

```

---

### False Position Method

#### False Position Theory
[Add your theory content here]

#### False Position Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double evalPoly(const vector<double>& coef, double x)
{
    double val = 0;
    int deg = coef.size() - 1;
    for(int i = 0; i <= deg; i++)
        val += coef[i] * pow(x, deg - i);
    return val;
}

string polyStr(const vector<double>& coef)
{
    stringstream ss;
    int deg = coef.size() - 1;
    bool first = true;

    for(int i = 0; i <= deg; i++)
    {
        double c = coef[i];
        int p = deg - i;
        if(fabs(c) < 1e-10) continue;

        if(!first) ss << (c > 0 ? " + " : " - ");
        else if(c < 0) ss << "-";

        double v = fabs(c);

        if(p == 0) ss << fixed << setprecision(4) << v;
        else if(p == 1) ss << fixed << setprecision(4) << v << "*x";
        else ss << fixed << setprecision(4) << v << "*x^" << p;

        first = false;
    }
    return ss.str();
}

double searchLimit(const vector<double>& coef)
{
    if(coef.size() < 3) return 10.0;
    double a = coef[0], b = coef[1], c = coef[2];
    return sqrt(abs((b / a) * (b / a) - 2 * (c / a)));
}

int main()
{
    ifstream fin("Input.txt");
    ofstream fout("Output.txt");

    if(!fin.is_open())
    {
        cerr << "Error: cannot open Input.txt\n";
        return 1;
    }

    int T;
    fin >> T;

    for(int t = 1; t <= T; t++)
    {
        int deg;
        fin >> deg;
        vector<double> coef(deg + 1);
        for(int i = 0; i <= deg; i++) fin >> coef[i];

        double tol;
        fin >> tol;

        cout << "TestCase " << t << "\n";
        fout << "TestCase " << t << "\n";

        string pstr = polyStr(coef);
        cout << "Polynomial: f(x)=" << pstr << "\n";
        fout << "Polynomial: f(x)=" << pstr << "\n";

        cout << "Tolerance: " << scientific << tol << "\n";
        fout << "Tolerance: " << scientific << tol << "\n";

        double lim = min(searchLimit(coef), 100.0);
        cout << "Search range: [" << -lim << ", " << lim << "]\n";
        fout << "Search range: [" << -lim << ", " << lim << "]\n";

        vector<double> roots;
        int rootCount = 0;

        for(double x = -lim; x < lim; x += 0.5)
        {
            double xl = x, xr = x + 0.5;
            double fl = evalPoly(coef, xl);
            double fr = evalPoly(coef, xr);

            if(fabs(fl) < tol) { roots.push_back(xl); continue; }
            if(fabs(fr) < tol) { roots.push_back(xr); continue; }

            if(fl * fr < 0)
            {
                int iter = 0;
                double xm;
                while(fabs(xr - xl) > tol)
                {
                    iter++;
                    xm = (xl * fr - xr * fl) / (fr - fl);
                    double fm = evalPoly(coef, xm);

                    if(fabs(fm) < tol) break;

                    if(fm * fl < 0) { xr = xm; fr = fm; }
                    else { xl = xm; fl = fm; }
                }

                roots.push_back(xm);
                rootCount++;
                cout << "Root " << rootCount << " in [" << fixed << setprecision(4) << xl << ", " << xr << "] = " << xm << " (iter: " << iter << ")\n";
                fout << "Root " << rootCount << " in [" << fixed << setprecision(4) << xl << ", " << xr << "] = " << xm << " (iter: " << iter << ")\n";
            }
        }

        sort(roots.begin(), roots.end());
        roots.erase(unique(roots.begin(), roots.end(), [tol](double a, double b){ return fabs(a-b) < tol; }), roots.end());

        cout << "All Real Roots:\n";
        fout << "All Real Roots:\n";

        if(roots.empty()) { cout << "None\n"; fout << "None\n"; }
        else
        {
            for(size_t i = 0; i < roots.size(); i++)
            {
                cout << "Root " << i+1 << ": " << fixed << setprecision(6) << roots[i] << "\n";
                fout << "Root " << i+1 << ": " << fixed << setprecision(6) << roots[i] << "\n";
            }
        }

        cout << "\n";
        fout << "\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### False Position Input
```
3
3
1 -6 11 -6
0.001
2
1 -4 4
0.001
5
1 -3 -7 27 -10 -8
0.001

```

#### False Position Output
```
TestCase 1
Polynomial: f(x)=1.0000*x^3-6.0000*x^2+11.0000*x-6.0000
Tolerance: 1.000000e-03
Search range: [-3.0, 3.0]
Root 1 in [0.5000, 1.0000] = 1.000000 (iter: 11)
Root 2 in [1.5000, 2.0000] = 2.000000 (iter: 12)
Root 3 in [2.5000, 3.0000] = 3.000000 (iter: 12)
All Real Roots:
Root 1: 1.000000
Root 2: 2.000000
Root 3: 3.000000

TestCase 2
Polynomial: f(x)=1.0000*x^2-4.0000*x+4.0000
Tolerance: 1.000000e-03
Search range: [-10.0, 10.0]
Root 1 in [1.5000, 2.0000] = 2.000000 (iter: 13)
All Real Roots:
Root 1: 2.000000

TestCase 3
Polynomial: f(x)=1.0000*x^5-3.0000*x^4-7.0000*x^3+27.0000*x^2-10.0000*x-8.0000
Tolerance: 1.000000e-03
Search range: [-5.47723, 5.47723]
Root 1 in [-2.0000, -1.5000] = -2.000244 (iter: 14)
Root 2 in [0.5000, 1.0000] = 0.999756 (iter: 15)
Root 3 in [1.5000, 2.0000] = 2.000244 (iter: 15)
All Real Roots:
Root 1: -2.000244
Root 2: 0.999756
Root 3: 2.000244

```

---

### Secant Method

#### Secant Theory
[Add your theory content here]

#### Secant Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> foundRoots;
ofstream fout;

double evalPoly(double x, const vector<double>& coef)
{
    double sum = 0;
    for(int i = 0; i < coef.size(); i++)
        sum += coef[i] * pow(x, i);
    return sum;
}

void printPolynomial(const vector<double>& coef)
{
    bool firstTerm = true;
    int n = coef.size();

    for(int i = n-1; i >= 0; i--)
    {
        if(fabs(coef[i]) < 1e-12) continue;

        if(!firstTerm)
        {
            if(coef[i] < 0) cout << " - ", fout << " - ";
            else cout << " + ", fout << " + ";
        }
        else
        {
            if(coef[i] < 0) cout << "-", fout << "-";
            firstTerm = false;
        }

        double c = fabs(coef[i]);

        if(i == 0) cout << c, fout << c;
        else if(i == 1) cout << c << "x", fout << c << "x";
        else cout << c << "x^" << i, fout << "x^" << i;
    }

    cout << endl;
    fout << endl;
}

void secantMethod(double x0, double x1, double tol, const vector<double>& coef)
{
    double xNext;
    int iterations = 0;

    while(fabs(x1 - x0) >= tol && fabs(evalPoly(x1, coef)) >= tol)
    {
        double f0 = evalPoly(x0, coef);
        double f1 = evalPoly(x1, coef);

        if(fabs(f1 - f0) < 1e-12) return;

        xNext = (x0 * f1 - x1 * f0) / (f1 - f0);
        iterations++;

        x0 = x1;
        x1 = xNext;
    }

    if(fabs(evalPoly(xNext, coef)) > tol) return;

    for(double r : foundRoots)
        if(fabs(r - xNext) < tol) return;

    foundRoots.push_back(xNext);

    cout << fixed << setprecision(6);
    fout << fixed << setprecision(6);

    cout << "Interval: [" << x0 << ", " << x1 << "] Root: " << xNext
         << " Iterations: " << iterations << endl;
    fout << "Interval: [" << x0 << ", " << x1 << "] Root: " << xNext
         << " Iterations: " << iterations << endl;
}

int main()
{
    ifstream fin("Input.txt");
    fout.open("Output.txt");

    if(!fin.is_open())
    {
        cout << "Cannot open Input.txt" << endl;
        return 0;
    }

    int testCases;
    fin >> testCases;

    while(testCases--)
    {
        foundRoots.clear();

        int deg;
        fin >> deg;
        vector<double> coef(deg + 1);
        for(double &x : coef) fin >> x;
        reverse(coef.begin(), coef.end());

        cout << "\nPolynomial: ";
        fout << "\nPolynomial: ";
        printPolynomial(coef);

        double xmax = 0;
        for(int i = 0; i <= deg; i++)
            xmax = max(xmax, 1.0 + fabs(coef[i]/coef[deg]));

        cout << "Root bound: " << xmax << endl;
        fout << "Root bound: " << xmax << endl;

        double step = 0.45;
        double tol = 1e-3;

        for(double i = -xmax; i <= xmax; i += step)
            secantMethod(i, i + step, tol, coef);
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Secant Input
```
2
3
1 -6 11 -6
5
1 -3 -7 27 -10 -8

```

#### Secant Output
```
Polynomial: 1x^3 - 6x^2 + 11x - 6
Root bound: 3
Interval: [0.5000, 1.0000] Root: 1.000000 Iterations: 6
Interval: [1.5000, 2.0000] Root: 2.000000 Iterations: 7
Interval: [2.5000, 3.0000] Root: 3.000000 Iterations: 7

Polynomial: 1x^5 - 3x^4 - 7x^3 + 27x^2 - 10x - 8
Root bound: 5.47723
Interval: [-2.0000, -1.5000] Root: -2.000244 Iterations: 8
Interval: [0.5000, 1.0000] Root: 0.999756 Iterations: 9
Interval: [1.5000, 2.0000] Root: 2.000244 Iterations: 9

```

---

### Newton Raphson Method

#### Newton Raphson Theory
[Add your theory content here]

#### Newton Raphson Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double polyVal(const vector<double>& coef, double x)
{
    double sum = 0;
    int n = coef.size();
    for(int i = 0; i < n; i++)
        sum += coef[i] * pow(x, n-1-i);
    return sum;
}

double polyDeriv(const vector<double>& coef, double x)
{
    double sum = 0;
    int n = coef.size();
    for(int i = 0; i < n-1; i++)
        sum += (n-1-i) * coef[i] * pow(x, n-2-i);
    return sum;
}


string polyToStr(const vector<double>& coef)
{
    stringstream ss;
    int n = coef.size();
    bool first = true;

    for(int i = 0; i < n; i++)
    {
        double v = coef[i];
        int p = n-1-i;
        if(fabs(v) < 1e-12) continue;

        if(!first && v >= 0) ss << " + ";
        else if(v < 0) { ss << " - "; v = -v; }

        if(p == 0) ss << v;
        else if(p == 1) ss << v << "x";
        else ss << v << "x^" << p;

        first = false;
    }
    return ss.str();
}

int main()
{
    ifstream fin("Input.txt");
    ofstream fout("Output.txt");

    if(!fin.is_open())
    {
        cerr << "Error: cannot open Input.txt\n";
        return 1;
    }

    int T;
    fin >> T;

    for(int tc = 1; tc <= T; tc++)
    {
        int deg;
        fin >> deg;
        vector<double> coef(deg+1);
        for(int i = 0; i <= deg; i++) fin >> coef[i];

        double start, end, tol;
        fin >> start >> end >> tol;

        int maxIter;
        fin >> maxIter;

        cout << "TestCase " << tc << "\n";
        fout << "TestCase " << tc << "\n";

        string pstr = polyToStr(coef);
        cout << "Polynomial: f(x)=" << pstr << "\n";
        fout << "Polynomial: f(x)=" << pstr << "\n";

        cout << "Search interval: [" << start << ", " << end << "]\n";
        fout << "Search interval: [" << start << ", " << end << "]\n";

        cout << "Tolerance: " << scientific << tol << "\n";
        fout << "Tolerance: " << scientific << tol << "\n";

        cout << "Max iterations: " << maxIter << "\n";
        fout << "Max iterations: " << maxIter << "\n";

        set<double> roots;
        int count = 1;

        for(double x = start; x <= end; x += 0.5)
        {
            double curr = x, next;
            int iter;

            for(iter = 1; iter <= maxIter; iter++)
            {
                double fx = polyVal(coef, curr);
                double fdx = polyDeriv(coef, curr);

                if(fabs(fdx) < 1e-12) break;

                next = curr - fx / fdx;

                if(fabs(next - curr) < tol || fabs(fx) < tol) break;

                curr = next;
            }

            double r = round(next * 1000.0) / 1000.0;
            if(!isfinite(r)) continue;

            if(roots.find(r) == roots.end())
            {
                roots.insert(r);
                cout << "Root " << count << " = " << fixed << setprecision(6) << r
                     << " (Iterations: " << iter << ")\n";
                fout << "Root " << count << " = " << fixed << setprecision(6) << r
                     << " (Iterations: " << iter << ")\n";
                count++;
            }
        }

        if(roots.empty())
        {
            cout << "No real roots found.\n";
            fout << "No real roots found.\n";
        }

        cout << "\n";
        fout << "\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Newton Raphson Input
```
2
3
1 -6 11 -6
0 3 0.001
50
2
1 -4 4
0 3 0.001
50

```

#### Newton Raphson Output
```
TestCase 1
Polynomial: f(x)=1x^3 - 6x^2 + 11x - 6
Search interval: [0, 3]
Tolerance: 1.000000e-03
Max iterations: 50
Root 1 = 1.000000 (Iterations: 5)
Root 2 = 2.000000 (Iterations: 5)
Root 3 = 3.000000 (Iterations: 5)

TestCase 2
Polynomial: f(x)=1x^2 - 4x + 4
Search interval: [0, 3]
Tolerance: 1.000000e-03
Max iterations: 50
Root 1 = 2.000000 (Iterations: 6)

```

---

### Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory
[Add your theory content here]

#### Gauss Elimination Code
```python
# Add your code here
```

#### Gauss Elimination Input
```
[Add your input format here]
```

#### Gauss Elimination Output
```
[Add your output format here]
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
[Add your theory content here]

#### Gauss Jordan Code
```python
# Add your code here
```

#### Gauss Jordan Input
```
[Add your input format here]
```

#### Gauss Jordan Output
```
[Add your output format here]
```

---

### LU Decomposition Method

#### LU Decomposition Theory
[Add your theory content here]

#### LU Decomposition Code
```python
# Add your code here
```

#### LU Decomposition Input
```
[Add your input format here]
```

#### LU Decomposition Output
```
[Add your output format here]
```

---

### Matrix Inversion

#### Matrix Inversion Theory
[Add your theory content here]

#### Matrix Inversion Code
```python
# Add your code here
```

#### Matrix Inversion Input
```
[Add your input format here]
```

#### Matrix Inversion Output
```
[Add your output format here]
```

---

### Solution of Differential Equations

### Runge-Kutta 4th Order Method

#### Runge-Kutta 4th Order Theory
[Add your theory content here]

#### Runge-Kutta 4th Order Code
```python
# Add your code here
```

#### Runge-Kutta 4th Order Input
```
[Add your input format here]
```

#### Runge-Kutta 4th Order Output
```
[Add your output format here]
```

---

### Numerical Integration

### Simpson's 1/3 Rule

#### Simpson's 1/3 Theory
[Add your theory content here]

#### Simpson's 1/3 Code
```cpp
#include <bits/stdc++.h>
using namespace std;

using db = double;

db evalPoly(db x, const vector<db>& coeff)
{
    db sum = 0;
    db powx = 1;
    for(db c : coeff)
    {
        sum += c * powx;
        powx *= x;
    }
    return sum;
}

void simpsonRule(int tc, ifstream &fin, ofstream &fout)
{
    int deg;
    fin >> deg;
    vector<db> coef(deg+1);
    for(int i=0; i<=deg; i++) fin >> coef[i];

    db a, b;
    int n;
    fin >> a >> b >> n;

    if(n % 2 != 0)
    {
        cout << "Test Case #" << tc << ": Error - n must be even!\n\n";
        fout << "Test Case #" << tc << ": Error - n must be even!\n\n";
        return;
    }

    db h = (b - a) / n;
    db sum = evalPoly(a, coef) + evalPoly(b, coef);

    for(int i=1; i<n; i++)
    {
        db x = a + i*h;
        sum += (i % 2 == 1) ? 4*evalPoly(x, coef) : 2*evalPoly(x, coef);
    }

    db result = (h/3) * sum;

    cout << "Test Case #" << tc << "\n";
    cout << "Degree: " << deg << "\n";
    cout << "Coefficients: ";
    for(db c : coef) cout << c << " ";
    cout << "\nInterval: [" << a << ", " << b << "]\n";
    cout << "Intervals (n): " << n << "\n";
    cout << "Step size (h): " << h << "\n";
    cout << "Integral result: " << result << "\n\n";

    fout << "Test Case #" << tc << "\n";
    fout << "Degree: " << deg << "\n";
    fout << "Coefficients: ";
    for(db c : coef) fout << c << " ";
    fout << "\nInterval: [" << a << ", " << b << "]\n";
    fout << "Intervals (n): " << n << "\n";
    fout << "Step size (h): " << h << "\n";
    fout << "Integral result: " << result << "\n\n";
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin)
    {
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    int t;
    fin >> t;
    cout << "Total Test Cases: " << t << "\n\n";
    fout << "Total Test Cases: " << t << "\n\n";

    for(int tc = 1; tc <= t; tc++)
        simpsonRule(tc, fin, fout);

    cout << "All results written to output.txt!\n";
    return 0;
}

```

#### Simpson's 1/3 Input
```
2
3
1 0 -2 1
0 2 4
2
2
2 -3 1
-1 1 6

```

#### Simpson's 1/3 Output
```
Test Case #1
Degree: 3
Coefficients: 1 0 -2 1 
Interval: [0, 2]
Intervals (n): 4
Step size (h): 0.5
Integral result: 0.666667

Test Case #2
Degree: 2
Coefficients: 2 -3 1 
Interval: [-1, 1]
Intervals (n): 6
Step size (h): 0.333333
Integral result: 0.000000

```

---

### Simpson's 3/8 Rule

#### Simpson's 3/8 Theory
[Add your theory content here]

#### Simpson's 3/8 Code
```python
#include <bits/stdc++.h>
using namespace std;

using db = double;

db evalPoly(db x, const vector<db>& coeff)
{
    db sum = 0, powx = 1;
    for(db c : coeff)
    {
        sum += c * powx;
        powx *= x;
    }
    return sum;
}

void simpson38(int tc, ifstream &fin, ofstream &fout)
{
    int deg;
    fin >> deg;
    vector<db> coef(deg + 1);
    for(int i=0; i<=deg; i++) fin >> coef[i];

    db a, b;
    int n;
    fin >> a >> b >> n;

    if(n % 3 != 0)
    {
        cout << "Test Case #" << tc << ": Error - n must be a multiple of 3!\n\n";
        fout << "Test Case #" << tc << ": Error - n must be a multiple of 3!\n\n";
        return;
    }

    db h = (b - a) / n;
    db sum = evalPoly(a, coef) + evalPoly(b, coef);

    for(int i=1; i<n; i++)
    {
        db x = a + i*h;
        sum += (i % 3 == 0) ? 2*evalPoly(x, coef) : 3*evalPoly(x, coef);
    }

    db result = (3*h/8) * sum;

    cout << "Test Case #" << tc << "\n";
    cout << "Degree: " << deg << "\n";
    cout << "Coefficients: ";
    for(db c : coef) cout << c << " ";
    cout << "\nInterval: [" << a << ", " << b << "]\n";
    cout << "Intervals (n): " << n << "\n";
    cout << "Step size (h): " << h << "\n";
    cout << "Integral result: " << result << "\n\n";

    fout << "Test Case #" << tc << "\n";
    fout << "Degree: " << deg << "\n";
    fout << "Coefficients: ";
    for(db c : coef) fout << c << " ";
    fout << "\nInterval: [" << a << ", " << b << "]\n";
    fout << "Intervals (n): " << n << "\n";
    fout << "Step size (h): " << h << "\n";
    fout << "Integral result: " << result << "\n\n";
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin)
    {
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    int t;
    fin >> t;
    cout << "Total Test Cases: " << t << "\n\n";
    fout << "Total Test Cases: " << t << "\n\n";

    for(int tc=1; tc<=t; tc++)
        simpson38(tc, fin, fout);

    cout << "All results written to output.txt!\n";
    return 0;
}

```

#### Simpson's 3/8 Input
```
2
3
1 0 -2 1
0 2 6
2
2
2 -3 1
-1 1 3

```

#### Simpson's 3/8 Output
```
Test Case #1
Degree: 3
Coefficients: 1 0 -2 1 
Interval: [0, 2]
Intervals (n): 6
Step size (h): 0.333333
Integral result: 0.666667

Test Case #2
Degree: 2
Coefficients: 2 -3 1 
Interval: [-1, 1]
Intervals (n): 3
Step size (h): 0.666667
Integral result: 0.000000

```

---


