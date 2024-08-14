// See https://aka.ms/new-console-template for more information

using static System.Math;

public class PreProcess
{


    public static List<float> EnsureDesiredLength(List<float> signal, int desiredLength)
    {
        int currentLength = signal.Count;
        //Console.WriteLine($"Current length of the signal: {currentLength}");

        if (currentLength < desiredLength)
        {
            float lastValue = signal[currentLength - 1];
            //Console.WriteLine($"Last value in the current signal: {lastValue}");

            for (int i = currentLength; i < desiredLength; i++)
            {
                signal.Add(lastValue);
            }

            //Console.WriteLine($"Signal length after filling: {signal.Count}");
        }
        return signal;
    }

    public static List<List<float>> SlidingWindow(List<float> ecgSignal, int windowSize)
    {
        List<List<float>> segments = new List<List<float>>();
        int start = 0;
        int end = windowSize;

        while (end <= ecgSignal.Count)
        {
            List<float> segment = ecgSignal.GetRange(start, windowSize);
            segments.Add(segment);
            start += windowSize ;
            end += windowSize ;
        }

        return segments;
    }
    public static List<float> MinMax(List<float> ecgSignal)
    {
        List<float> MinMaxedSignal = new List<float>();
        float min = ecgSignal.Min();
        float max = ecgSignal.Max();
        
        for (int i = 0; i < ecgSignal.Count; i++)
        {
            float normalizedValue = 2 * ((ecgSignal[i] - min) / (max - min)) - 1;
            MinMaxedSignal.Add(normalizedValue);
        }
        
        return MinMaxedSignal;
    }

    public static double CubicHermitePoly(double t, double startPoint, double endPoint, double p0, double p1, double m0, double m1)
    {
        double t2 = t * t;
        double t3 = t2 * t;
        return (2*t3 - 3*t2 + 1)*p0 + (t3 - 2*t2 + t)*m0*(endPoint-startPoint) + (-2*t3 + 3*t2)*p1 + (t3 - t2)*m1*(endPoint-startPoint); 
    }

    public static double[] GetTangents(double fs, double[] Y)
    {
        int outLength = Y.Length;
        double[] outM = new double[Y.Length];
        
        outM[0] = fs * (Y[1] - Y[0]);
        outM[outLength - 1] = fs * (Y[outLength - 1] - Y[outLength - 2]);

        for (int i = 1; i < outLength - 1; i++)
        {
            outM[i] = 0.5 * fs * (Y[i + 1] - Y[i-1]);
        }
        
        return outM;
    }

    public static double[] CubicHermiteInterpol(double[] Y, double fsOrig, double fsTarget) // resampling function
    {
        double ratio = fsTarget / fsOrig;
        int outLength = (int)((Y.Length - 1) * ratio);
        double[] outVector = new double[outLength];
        double[] Tangents = GetTangents(fsOrig, Y);
        for (int i = 0; i < outLength; i++) 
        {

           

            int LowSample = (int)(i / ratio);
            int Highsample = LowSample + 1;

            
                   
                    double p0 = Y[LowSample];
                    double p1 = Y[Highsample];
                    double m0 = Tangents[LowSample];
                    double m1 = Tangents[Highsample];
                    outVector[i] = (CubicHermitePoly(i / fsTarget - LowSample / fsOrig, LowSample / fsOrig, Highsample / fsOrig, p0, p1, m0, m1));
                
            

        }

        return outVector;

    }
}

