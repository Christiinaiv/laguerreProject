using System;
using System.Collections.Generic;
using System.Globalization;
using static System.Net.Mime.MediaTypeNames;


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
            if (n == 0) return Math.Sqrt(Sigma) * Math.Exp((-Beta*t)/2);
            else if (n == 1) return Math.Sqrt(Sigma) * (1-Sigma*t) * Math.Exp((-Beta*t)/2);
            else
            {
                var l0 = Math.Sqrt(Sigma) * Math.Exp((-Beta * t) / 2);
                var l1 = Math.Sqrt(Sigma) * (1 - Sigma * t) * Math.Exp((-Beta * t) / 2);
                double li = 0;
                for (var i = 2; i <= n; i++)
                {
                    li = ((2 * i - 1 - Sigma * t) / (double)i) * l1 - ((i - 1) / (double)i) * l0;
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

    public IntegrateLaguerre(Laguerre laguerre, int maxSteps = 10000000)
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
            double deltaT = T / steps;

            var tValues = new List<double>();
            for (int i = 0; i <= steps; i++)
            {
                tValues.Add(i * deltaT);
            }

            var lkt = new List<double>();
            var ft = new List<double>();
            for (int i = 0; i < tValues.Count; i++)
            {
                double t = tValues[i];
                lkt.Add(ComputeLaguerre(k, t));
                ft.Add(f(t));
            }

            double integral = 0;
            for (int i = 0; i < tValues.Count; i++)
            {
                integral += ft[i] * lkt[i] * deltaT * Math.Exp(-Alpha * tValues[i]);
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

 
class Prog
{
    static Func<double, double> ParseFunction(string expression)
    {
        var tokens = Tokenize(expression);
        var rpn = ToRPN(tokens);
        return t => EvaluateRPN(rpn, t);
    }

    static List<string> Tokenize(string expr)
    {
        var tokens = new List<string>();
        int i = 0;
        while (i < expr.Length)
        {
            char c = expr[i];

            if (char.IsWhiteSpace(c))
            {
                i++;
                continue;
            }

            if (char.IsDigit(c) || c == '.')
            {
                int start = i;
                while (i < expr.Length && (char.IsDigit(expr[i]) || expr[i] == '.')) i++;

                if (i < expr.Length && expr[i] == '/')
                {
                    i++;
                    int startDenom = i;
                    while (i < expr.Length && (char.IsDigit(expr[i]) || expr[i] == '.')) i++;
                    string numerator = expr.Substring(start, startDenom - start - 1);
                    string denominator = expr.Substring(startDenom, i - startDenom);
                    double num = double.Parse(numerator, CultureInfo.InvariantCulture);
                    double den = double.Parse(denominator, CultureInfo.InvariantCulture);
                    tokens.Add((num / den).ToString(CultureInfo.InvariantCulture));
                }
                else
                {
                    tokens.Add(expr.Substring(start, i - start));
                }
            }
            else if (char.IsLetter(c))
            {
                int start = i;
                while (i < expr.Length && char.IsLetter(expr[i])) i++;
                tokens.Add(expr.Substring(start, i - start));
            }
            else if (c == '-' && (tokens.Count == 0 || "+-*/^(".Contains(tokens[^1])))
            {
                tokens.Add("u-"); // унарний мінус
                i++;
            }
            else if ("+-*/^(),".Contains(c))
            {
                tokens.Add(c.ToString());
                i++;
            }
            else
            {
                throw new Exception($"Невідомий символ: {c}");
            }
        }
        return tokens;
    }

    static List<string> ToRPN(List<string> tokens)
    {
        var output = new List<string>();
        var opStack = new Stack<string>();
        var precedence = new Dictionary<string, int>
        {
            ["u-"] = 5,
            ["^"] = 4,
            ["*"] = 3,
            ["/"] = 3,
            ["+"] = 2,
            ["-"] = 2
        };

        foreach (var token in tokens)
        {
            if (double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out _) || token == "t")
            {
                output.Add(token);
            }
            else if (IsFunction(token))
            {
                opStack.Push(token);
            }
            else if (token == ",")
            {
                while (opStack.Peek() != "(")
                    output.Add(opStack.Pop());
            }
            else if (precedence.ContainsKey(token))
            {
                while (opStack.Count > 0 && precedence.ContainsKey(opStack.Peek()) && precedence[opStack.Peek()] >= precedence[token])
                {
                    output.Add(opStack.Pop());
                }
                opStack.Push(token);
            }
            else if (token == "(")
            {
                opStack.Push(token);
            }
            else if (token == ")")
            {
                while (opStack.Count > 0 && opStack.Peek() != "(")
                {
                    output.Add(opStack.Pop());
                }
                if (opStack.Count == 0 || opStack.Peek() != "(")
                    throw new Exception("Дужки не збігаються");
                opStack.Pop();

                if (opStack.Count > 0 && IsFunction(opStack.Peek()))
                    output.Add(opStack.Pop());
            }
            else
            {
                throw new Exception($"Невідомий токен: {token}");
            }
        }

        while (opStack.Count > 0)
        {
            if (opStack.Peek() == "(")
                throw new Exception("Дужки не збігаються");
            output.Add(opStack.Pop());
        }

        return output;
    }

    static double EvaluateRPN(List<string> rpn, double t)
    {
        var stack = new Stack<double>();
        foreach (var token in rpn)
        {
            if (double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out double val))
            {
                stack.Push(val);
            }
            else if (token == "t")
            {
                stack.Push(t);
            }
            else if (token == "u-")
            {
                if (stack.Count < 1) throw new Exception("Недостатньо аргументів для унарного мінуса.");
                double x = stack.Pop();
                stack.Push(-x);
            }
            else if (IsOperator(token))
            {
                if (stack.Count < 2) throw new Exception("Недостатньо аргументів.");
                double b = stack.Pop();
                double a = stack.Pop();

                stack.Push(token switch
                {
                    "+" => a + b,
                    "-" => a - b,
                    "*" => a * b,
                    "/" => a / b,
                    "^" => Math.Pow(a, b),
                    _ => throw new Exception($"Невідомий оператор: {token}")
                });
            }
            else if (IsFunction(token))
            {
                if (stack.Count < 1) throw new Exception("Недостатньо аргументів для функції.");
                double x = stack.Pop();
                stack.Push(token switch
                {
                    "sin" => Math.Sin(x),
                    "cos" => Math.Cos(x),
                    "tan" => Math.Tan(x),
                    "exp" => Math.Exp(x),
                    _ => throw new Exception($"Невідома функція: {token}")
                });
            }
            else
            {
                throw new Exception($"Невідомий токен при обчисленні: {token}");
            }
        }

        if (stack.Count != 1)
            throw new Exception("Некоректний вираз.");

        return stack.Pop();
    }

    static bool IsOperator(string token) => token is "+" or "-" or "*" or "/" or "^";
    static bool IsFunction(string token) => token is "sin" or "cos" or "tan" or "exp";



    static void Main()
        {
            var ifilePath = ".\\input_params.csv";
            var ofilePath = ".\\output_results.csv";
            using (StreamReader reader = new StreamReader(ifilePath))
            {
                reader.ReadLine(); 

                string line = reader.ReadLine();
                if (string.IsNullOrEmpty(line))
                {
                    throw new InvalidOperationException("Файл input_params.csv порожній або не містить даних.");
                }

                string[] values = line.Split(',');
            
                if (values.Length < 8)
                {
                    throw new InvalidOperationException($"Некоректна кількість стовпців у файлі input_params.csv. Очікується 8, отримано {values.Length}.");
                }

                var t_values = values[5].Split(';');
            
                if (!int.TryParse(values[7], out int k))
                {
                    throw new FormatException($"Неможливо перетворити значення k ({values[7]}) у ціле число.");
                }

                // Зчитуємо функцію з файлу
                string function = values[6];
                Func<double, double> f;
                try
                {
                    f = ParseFunction(function);

            }
                catch (ArgumentException ex)
                {
                    throw new InvalidOperationException($"Помилка при парсингу функції: {ex.Message}");
                }

                var parsedT = double.Parse(values[0], CultureInfo.InvariantCulture);
                var parsedN = int.Parse(values[1]);
                var parsedBeta = double.Parse(values[2], CultureInfo.InvariantCulture);
                var parsedSigma = double.Parse(values[3], CultureInfo.InvariantCulture);
                var parsedEpsilon = double.Parse(values[4], CultureInfo.InvariantCulture);
                var parsedTValues = new List<double>();
                foreach (var tVal in t_values)
                {
                    if (double.TryParse(tVal, NumberStyles.Float, CultureInfo.InvariantCulture, out double parsedValue))
                    {
                        parsedTValues.Add(parsedValue);
                    }
                    else
                    {
                        throw new FormatException($"Некоректний формат значення t: {tVal}");
                    }
                }

                var test1 = new Laguerre(parsedT, parsedN, parsedBeta, parsedSigma, parsedEpsilon, parsedTValues);
                var test2 = new IntegrateLaguerre(test1);
                var test3 = new ReverseLaguerre(test1, new List<double>() { });

                var table = test1.TabulateLaguerre();
                using (StreamWriter writer = new StreamWriter(ofilePath))
                {
                    writer.WriteLine("integral_result,coefficients,reverse_value,t_epsilon");
                    writer.WriteLine(
                        $"{test2.Integrate(f, k).ToString(CultureInfo.InvariantCulture)}," + 
                        $"\"{string.Join(";", test3.LaguerreCoefficients(test2, f).Select(c => c.ToString(CultureInfo.InvariantCulture)))}\"," +
                        $"{test3.ReverseFunc(parsedT).ToString(CultureInfo.InvariantCulture)}," +
                        $"{test1.FindTEpsilon(k).ToString(CultureInfo.InvariantCulture)}"
                    );
                }

                System.Console.WriteLine("ALL_good");
            }
        }
}