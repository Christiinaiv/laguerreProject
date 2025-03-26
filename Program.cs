using System;
using System.Collections.Generic;


class Laguerre
{
    public double T { get; set; }
    public int N { get; set; }
    public double Beta { get; set; }
    public double Sigma { get; set; }
    public double Epsilon { get; set; }
    public List<double> TValues { get; set; }
    protected double Alpha => Sigma - Beta;

    public Laguerre(double t, int n, double beta, double sigma, double epsilon, List<double> tValues)
    {
        T = t;
        N = n;
        Beta = beta;
        Sigma = sigma;
        Epsilon = epsilon;
        TValues = tValues;
    }


    public double ComputeLaguerre(int n, double t)
    {
        if (Beta >= 0 && Beta <= Sigma)
        {
            if (n == 0)
            {
                return Math.Sqrt(Sigma) * Math.Exp((-Beta * t) / 2);
            }
            else if (n == 1)
            {
                return Math.Sqrt(Sigma) * (1 - Sigma * t) * Math.Exp((-Beta * t) / 2);
            }
            else
            {
                var l0 = Math.Sqrt(Sigma) * Math.Exp((-Beta * t) / 2);
                var l1 = Math.Sqrt(Sigma) * (1 - Sigma * t) * Math.Exp((-Beta * t) / 2);
                double li = 0;
                for (var i = 2; i <= n; i++)
                {
                    li = ((2 * i - 1 - Sigma * t) / i) * l1 - ((i - 1) / (double)i) * l0;
                    l0 = l1;
                    l1 = li;
                }
                return li;
            }
        }
        else
        {
            throw new ArgumentException("Wrong data input! Beta should be in range [0;sigma]!");
        }
    }

    public Dictionary<string, List<double>> TabulateLaguerre()
    {
        var tabulate = new Dictionary<string, List<double>>();
        tabulate["t_val"] = TValues;

        for (var i = 0; i < N; i++)
        {
            var values = new List<double>();
            foreach (var t in TValues)
            {
                values.Add(ComputeLaguerre(i, t));
            }
            tabulate[$"L_{i}"] = values;
        }

        return tabulate;
    }


    public double FindTEpsilon(int n, double epsilon = 0.001, int maxT = 100, int steps = 1000)
    {
        var stepSize = (double)maxT / steps;
        for (double t = 0; t < maxT; t += stepSize)
        {
            bool allLess = true;

            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(ComputeLaguerre(i, t)) >= epsilon)
                {
                    allLess = false;
                    break;
                }
            }
            if (allLess)
            {
                return t;
            }
        }
        throw new ArgumentException("No value t");
    }
}

class IntegrateLaguerre : Laguerre
{
    public int MaxSteps { get; set; }

    public IntegrateLaguerre(Laguerre laguerre, int maxSteps = 10000)
        : base(laguerre.T, laguerre.N, laguerre.Beta, laguerre.Sigma, laguerre.Epsilon, laguerre.TValues)
    {
        MaxSteps = maxSteps;
    }

    public double Integrate(Func<double, double> f, int k)
    {
        int steps = 10;
        double prevIntegral = 0;
        double delta = double.PositiveInfinity;

        while (delta > Epsilon)
        {
            var tValues = new List<double>();
            double deltaT = T / steps;

            for (int i = 0; i <= steps; i++)
            {
                tValues.Add(i * deltaT);
            }

            var lkt = new List<double>();
            foreach (var t in tValues)
            {
                lkt.Add(ComputeLaguerre(k, t));
            }

            var integrand = new List<double>();
            for (int i = 0; i < tValues.Count; i++)
            {
                double ft = f(tValues[i]);
                integrand.Add(ft * lkt[i] * Math.Exp(-Alpha * tValues[i]));
            }

            double integral = 0;
            for (int i = 0; i < steps; i++)
            {
                integral += (integrand[i] + integrand[i + 1]) * deltaT / 2;
            }

            delta = Math.Abs(integral - prevIntegral);
            prevIntegral = integral;

            steps *= 2;
            if (steps > MaxSteps)
            {
                throw new InvalidOperationException("Maximum number of steps exceeded. Integral does not converge.");
            }
        }

        return prevIntegral;
    }
}

class ReverseLaguerre : Laguerre
{
    public List<double> Coefficients { get; set; }

    public ReverseLaguerre(Laguerre laguerre, List<double> coefficients = null)
        : base(laguerre.T, laguerre.N, laguerre.Beta, laguerre.Sigma, laguerre.Epsilon, laguerre.TValues)
    {
        if (coefficients != null)
        {
            Coefficients = coefficients;
        }
        else
        {
            Coefficients = new List<double>();
        }
    }

    public List<double> LaguerreCoefficients(IntegrateLaguerre integrator, Func<double, double> f)
    {
        var coeffs = new List<double>();
        for (int k = 0; k <= N; k++)
        {
            coeffs.Add(integrator.Integrate(f, k));
        }
        Coefficients = coeffs;
        return coeffs;
    }

    public double ReverseFunc(double t)
    {
        if (Coefficients.Count == 0)
        {
            throw new InvalidOperationException("Coefficients are not computed");
        }

        double result = 0;
        for (int k = 0; k < Coefficients.Count; k++)
        {
            result += Coefficients[k] * ComputeLaguerre(k, t);
        }
        return result;
    }

    public List<double> ConditionalFunc(List<double> tValues = null)
    {
        if (tValues == null)
        {
            tValues = new List<double>();
            double stepSize = T / 100;
            for (int i = 0; i <= 100; i++)
            {
                tValues.Add(i * stepSize);
            }
        }

        var results = new List<double>();
        foreach (var t in tValues)
        {
            results.Add(ReverseFunc(t));
        }
        return results;
    }
}


class Program
{
    static void Main()
    {
        var test = new Laguerre(5, 5, 2, 4, 0.001, new List<double>() { 0, 1, 2, 3, 4, 5 });
        var table = test.TabulateLaguerre();

        Console.Write("t_val\t");
        for (int i = 0; i < test.N; i++)
        {
            Console.Write($"L_{i}\t");
        }
        Console.WriteLine();

        for (int j = 0; j < table["t_val"].Count; j++)
        {
            Console.Write($"{table["t_val"][j]:0.###}\t");
            for (int i = 0; i < test.N; i++)
            {
                Console.Write($"{table[$"L_{i}"][j]:0.###}\t");
            }
            Console.WriteLine();
        }
        Console.WriteLine();
        Console.WriteLine($"T Epsilon: {test.FindTEpsilon(3)}");

        Console.WriteLine("-----------------------------");
        Console.WriteLine("Iнтегрування");
        var test2 = new IntegrateLaguerre(test);
        Console.WriteLine($"Iнтеграл: {test2.Integrate(x => Math.Sin(x), 3)}");

        var test3 = new ReverseLaguerre(test, new List<double>() { });
        Console.WriteLine($"Коефіцієнти: {string.Join(", ", test3.LaguerreCoefficients(test2, x => Math.Sin(x)))}");
        Console.WriteLine($"Реверсивне значення: {test3.ReverseFunc(5)}");
    }
}
