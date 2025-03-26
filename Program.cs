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
        if(Beta>=0 && Beta<=Sigma)
        {
            if(n==0)
            {
                return Math.Sqrt(Sigma) * Math.Exp((-Beta * t) / 2);
            }
            else if(n==1)
            {
                return Math.Sqrt(Sigma) * (1 - Sigma * t) * Math.Exp((-Beta * t) / 2);
            }
            else
            {
                var l0 = Math.Sqrt(Sigma) * Math.Exp((-Beta * t) / 2);
                var l1 = Math.Sqrt(Sigma) * (1 - Sigma * t) * Math.Exp((-Beta * t) / 2);
                double li = 0;
                for(var i=2; i<=n; i++)
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

        for(var i=0; i<N; i++)
        {
            var values = new List<double>();
            foreach(var t in TValues)
            {
                values.Add(ComputeLaguerre(i, t));
            }
            tabulate[$"L_{i}"] = values;
        }

        return tabulate;
    }


   public double FindTEpsilon(int n, double epsilon = 0.001, int maxT=100, int steps=1000)
    {
        var stepSize = maxT / steps;
        for (double t = 0; t < maxT; t += stepSize)
        {
            bool allLess = true;

            for (int i = 0; i < n; i++)
            {
                if(Math.Abs(ComputeLaguerre(i, t))>=epsilon)
                {
                    allLess = false;
                    break;
                }
            }
            if(allLess)
            {
                return t;
            }
        }
        throw new ArgumentException("No value t");
    }
}
