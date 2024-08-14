using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
namespace FilterTenet;

    public class AdaptiveMedianFilter<T> where T : IComparable<T>, IConvertible
    {
        private const int DEFAULT_HR = 60; // Replace with the appropriate default value.
        public const int FAKE_ECG_SAMPLE_VALUE = 40;
        T boundThr;
        T zeroT;
        T adcZeroT;
        T fakeT;
        public AdaptiveMedianFilter()
        {
            boundThr = (T)Convert.ChangeType(50, typeof(T));
            zeroT = (T)Convert.ChangeType(0, typeof(T));
            adcZeroT = (T)Convert.ChangeType(2048, typeof(T));
            fakeT = (T)Convert.ChangeType(FAKE_ECG_SAMPLE_VALUE, typeof(T));
        }

        public T[] Apply(T[] ecg, int fs, int[] rpeaks)
        {
            try
            {
                //Console.WriteLine("**************** Run AdaptiveMedianFilter ****************");
                //Console.WriteLine("input.length: " + ecg.Length + ", rpeaks.length: " + ((rpeaks != null && rpeaks.Length > 0) ? rpeaks.Length : 0));
                //Console.WriteLine("rpeaks: " + ((rpeaks != null && rpeaks.Length > 0) ? "[" + string.Join(",", rpeaks) + "]" : "--"));
                //Console.WriteLine("input: [" + string.Join(",", ecg) + "]");

                List<T> finalResult = new List<T>();

                const int windowDuration = 8; // default ECG segment size
                int windowSampleCount = fs * windowDuration;
                for (int i = 0; i < ecg.Length; i += windowSampleCount)
                {
                    int iTo = Math.Min(i + windowSampleCount, ecg.Length);
                    //Console.WriteLine("ecg.len: " + ecg.Length + ", i: " + i + ", iTo: " + iTo);

                    int medianWindowSize = GetMedianWindowSize(fs, rpeaks, i, iTo);

                    int compensateWidth = (medianWindowSize / 2) + 1;
                    int iWithCompensateLeft = Math.Max(0, i - compensateWidth);
                    int iToWithCompensateRight = Math.Min(ecg.Length, iTo + compensateWidth);
                    T[] ecgSegment = ecg.Skip(iWithCompensateLeft).Take(iToWithCompensateRight - iWithCompensateLeft).ToArray();

                    //Console.WriteLine("ecgSegment: [" + string.Join(",", ecgSegment) + "]");

                    // Run median filter for each peak using and append to final medianECg
                    // Get median array and subtract from the original signal
                    //Console.WriteLine("filterSignal2 call start");
                    List<T> medianResult = FilterSignal2(ecgSegment, medianWindowSize);
                    //Console.WriteLine("filterSignal2 call end");

                    // Return the result

                    int segmentPosfrom = 0;
                    int segmentPosTo = ecgSegment.Length;
                    if (iTo - i < ecgSegment.Length)
                    {
                        if (iWithCompensateLeft < i)
                        {
                            segmentPosfrom = (i - iWithCompensateLeft);
                        }
                        if (iTo < iToWithCompensateRight)
                        {
                            segmentPosTo = ecgSegment.Length - (iToWithCompensateRight - iTo);
                        }
                    }

                    //Console.WriteLine("------------------------------ecg.len: " + ecg.Length + "---------------------------------------");
                    //Console.WriteLine("i: " + i + ", iTo: " + iTo + ", iiToLen: " + (iTo - i) + ", ecgsegmentLen: " + ecgSegment.Length);
                    //Console.WriteLine("medianWindowSize: " + medianWindowSize + ", compensateWidth: " + compensateWidth);
                    //Console.WriteLine("iWithCompensateLeft: " + iWithCompensateLeft + ", iToWithCompensateRight: " + iToWithCompensateRight);
                    //Console.WriteLine("segmentPosfrom: " + segmentPosfrom + ", segmentPosTo: " + segmentPosTo + ", len: " + (segmentPosTo - segmentPosfrom));
                    //Console.WriteLine("-----------------------------------------------------------------------------------------");

                    finalResult.AddRange(medianResult.Skip(segmentPosfrom).Take(segmentPosTo - segmentPosfrom));
                }

                //Console.WriteLine("finalResultSize:" + finalResult.Count);
                //Console.WriteLine("output: [" + string.Join(",", finalResult) + "]");


                T avg = Utils.AvgIgnoreFake(ecg);

                T[] finalResultOnBaseline = new T[finalResult.Count];
                for (int i = 0; i < finalResult.Count; i++)
                {
                    bool isFake = fakeT.CompareTo(ecg[i]) == 0;
                    finalResultOnBaseline[i] = isFake
                        ? fakeT
                        : (T)(object)(Convert.ToSingle(avg) + Convert.ToSingle(finalResult[i]));
                }

                //Console.WriteLine("avg: " + avg + ", finalResult: [" + string.Join(",", finalResultOnBaseline) + "]");
                //WriteToCsv(finalResultOnBaseline, outputFilePath);
                return finalResultOnBaseline;
            }
            catch (Exception e)
            {
                Console.WriteLine("Crash: " + e.Message);
                return ecg;
            }
        }





        private int GetMedianWindowSize(int fs, int[] rpeaks, int rpeakFrom, int rpeakTo)
        {
            int sampleIntervalMsec = 1000 / fs;

            // Calc HR based on rPeaks on ECG segments
            int hrOnSegment = DEFAULT_HR;
            if (rpeaks != null && rpeaks.Length > 0)
            {
                try
                {
                    // Extract rpeaks of segments
                    List<int> rpeaksOnSegment = new List<int>();
                    for (int j = 0; j < rpeaks.Length; j++)
                    {
                        if (rpeakFrom <= rpeaks[j] && rpeaks[j] < rpeakTo)
                        {
                            rpeaksOnSegment.Add(rpeaks[j]);
                        }
                    }

                    // Find HR for each rpeak
                    List<int> hrOnRPeaksOnSegment = new List<int>();
                    for (int j = 1; j < rpeaksOnSegment.Count; j++)
                    {
                        hrOnRPeaksOnSegment.Add(GetMomentaryHR(sampleIntervalMsec, rpeaksOnSegment[j], rpeaksOnSegment[j - 1]));
                    }

                    // Find HRmean on segment
                    if(hrOnRPeaksOnSegment.Count()>0)
                        hrOnSegment = (int)hrOnRPeaksOnSegment.Average();

                }
                catch (Exception e)
                {
                    Console.WriteLine("Crash hrOnSegment: " + e.Message);
                }
            }

            if (hrOnSegment == 0)
            {
                Console.WriteLine("hrOnSegment is ZERO");
                hrOnSegment = DEFAULT_HR;
            }

            // Calculate median window size for each segment based on HR
            float flexibleWindowDur = ((60f / hrOnSegment) * 0.9f);
            int medianWindowSize = GetMedianFilterWidth(fs, flexibleWindowDur);
            //if (hrOnSegment == DEFAULT_HR)
                //Console.WriteLine("medianWindowSize: " + medianWindowSize + ", flexibleWindowDur: " + flexibleWindowDur);
            return medianWindowSize;
        }

        private int GetMomentaryHR(int sampleIntervalMsec, int rpeakCur, int rpeakPast)
        {
            return (int)(60000 / ((rpeakCur - rpeakPast) * sampleIntervalMsec));
        }

        private int GetMedianFilterWidth(int fs, float duration)
        {
            int sampleCount = (int)(fs * duration);
            sampleCount += ((sampleCount % 2) - 1);
            return sampleCount;
        }

        private T[] MedFilt2(T[] data, int window)
        {
            if (window % 2 != 1)
                throw new ArgumentException("Window size must be odd.");

            if (data == null || data.Length == 0)
                throw new ArgumentException("Data cannot be null or empty.");

            List<T> xList = data.ToList();
            T xListFirst = xList[0];
            T xListLast = xList[xList.Count - 1];

            List<T> result = new List<T>();

            int win2 = (int)((window - 1) / 2);
            for (int i = 0; i < data.Length; i++)
            {
                int from = Math.Max(0, i - win2);
                int to = Math.Min(i + win2 + 1, data.Length);
                List<T> sub = new List<T>(xList.GetRange(from, to - from));

                while (sub.Count < window)
                {
                    if (from == 0 && to == data.Length)
                    {
                        Console.WriteLine("subarray is too small!");
                    }
                    if (from != 0 && to != data.Length)
                    {
                        Console.WriteLine("subarray is suspicious!");
                    }

                    if (from == 0)
                    {
                        sub.Insert(0, xListFirst);
                    }
                    if (to == data.Length)
                    {
                        sub.Add(xListLast);
                    }
                }

                if (sub.Count != window)
                {
                    Console.WriteLine("subarray size is wrong");
                }
                T[] subArray = sub.OrderBy(item => item).ToArray();

                T foundMedianCast = typeof(T) == typeof(int)
                    ? (T)(object)Convert.ToInt32(subArray[win2])
                    : (T)(object)Convert.ToSingle(subArray[win2]);

                result.Add(foundMedianCast);
            }

            for (int p = 0; p < win2; p++)
            {
                int i = p;
                int j = result.Count - 1 - p;
                result[i] = result[win2];
                result[j] = result[result.Count - win2];
            }

            try
            {
                for (int p = 0; p < result.Count; p++)
                {
                    if (Convert.ToSingle(result[p]) == FAKE_ECG_SAMPLE_VALUE && Convert.ToSingle(data[p]) != FAKE_ECG_SAMPLE_VALUE)
                    {
                        bool found = false;
                        int pp = p - 1;
                        int pn = p + 1;
                        while (!found && (pp >= 0 || pn < result.Count))
                        {
                            bool foundInPrevious = pp >= 0 && Convert.ToSingle(result[pp]) != FAKE_ECG_SAMPLE_VALUE;
                            bool foundInNext = pn < result.Count && Convert.ToSingle(result[pn]) != FAKE_ECG_SAMPLE_VALUE;
                            if (foundInPrevious)
                            {
                                found = true;
                                result[p] = result[pp];
                            }
                            else if (foundInNext)
                            {
                                found = true;
                                result[p] = result[pn];
                            }
                            pp -= 1;
                            pn += 1;
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("Crash hrOnSegment22: " + e.Message);
            }

            return result.ToArray();
        }

        // get median array and subtract from org signal
        private List<T> FilterSignal2(T[] X, int m)
        {
            T[] X0 = new T[X.Length];
            for (int j = 0; j < X.Length; j++)
            {
                X0[j] = X[j];
            }
            X0 = MedFilt2(X0, m);
            List<T> X1 = new List<T>();

            for (int i = 0; i < X.Length; i++)
            {
                float sub = Convert.ToSingle(X[i]) - Convert.ToSingle(X0[i]);
                T subCast = typeof(T) == typeof(int)
                    ? (T)(object)Convert.ToInt32(sub)
                    : (T)(object)Convert.ToSingle(sub);
                X1.Add(subCast);
            }

            return X1;
        }

    }



    public class AdaptiveMovingAverage
    {
        private const string TAG = "AdaptiveMovingAverage";

        readonly float Af = 0.9F; // higher (0 < A_fall < 1) for slower cut-off fall
        readonly float Ar = 0.5F; // lower (0 < A_rise < 1) for faster cut-off rise
        static readonly float p =(float)Math.Sqrt(2); // Critically Damped: 2, Bessel: 3, Butterworth: sqrt(2)
        static readonly float g = 1f; // Critically Damped: 1, Bessel: 3, Butterworth: 1
       readonly float c = (float)Math.Sqrt(2f / (2 * g - p * p + Math.Sqrt(Math.Pow(2 * g - p * p, 2) + 4 * g * g)));


        float f0min = 4; // 2.5;  // min cut - off frequency
        static readonly float f0max = 40;   //no filter cut-off frequency
        float NFP = 0.8F;
        float NRP = 0.1F;
        int delmax = 6;


        readonly float pi = 3.1415926F;

        float f0 = f0max;
        float NM = 2.0F;

        public AdaptiveMovingAverage(float minCutOff = 4)
        {
            f0min = minCutOff;
            f0 = f0max;
            NM = 2.0F;
        }

        public float[] Apply_filter(float[] ecgArray, int fs)
        {           

            float[] previousEcgFAMA = null;
            float[] previousEcg = null;

            //Console.WriteLine(TAG + "which f0min: " + f0min);
            try
            {
                if (ecgArray == null || ecgArray.Length < 3)
                    return ecgArray;

                bool hasEnoughDataFromPrevious = previousEcg != null && previousEcg.Length > 2 && previousEcgFAMA != null && previousEcgFAMA.Length > 2;

                int Lecg = ecgArray.Length;
                float[] ecgF = new float[Lecg];
                ecgArray.CopyTo(ecgF, 0);

                for (int i = (hasEnoughDataFromPrevious ? 0 : 2); i < Lecg; i++)
                {
                    float fc = 0;

                    for (int j = 0; j <= delmax; j++)
                    {
                        fc = fc + 120 * f0min * Math.Abs(ecgArray[Math.Max(0, Math.Min(i + j, Lecg - 1))] - ((i - 1) < 0 && hasEnoughDataFromPrevious ? previousEcgFAMA[previousEcgFAMA.Length + (i - 1)] : ecgF[i - 1])) / (float)(delmax + 1);
                    }

                    fc = (float)Math.Max(Math.Max(f0min, f0 * Af), fc / NM);
                    f0 = (float)Math.Min(Math.Min(f0max, f0 / Ar), fc);

                    NM = (float)Math.Max(0.3, NM * (1 + Math.Min(NRP, (f0 / f0min) - 1 - NFP) / 100));

                    float w0 = (float)Math.Tan(pi * c * f0 / fs); // digital cut-off frequency
                    float a0 = g * w0 * w0 / (1 + p * w0 + g * w0 * w0);
                    float a1 = 2 * a0;
                    float a2 = a0;
                    float b1 = 2 * a0 * (1 / (g * w0 * w0) - 1);
                    float b2 = 1 - a0 - a1 - a2 - b1;

                    
                    int del = (int)Math.Min(Math.Min(delmax, Math.Round(0.24 * fs / f0)), Lecg - 1 - i);
                    try
                    {
                        ecgF[i] = a0 * ((i < 2 && hasEnoughDataFromPrevious && (i + del) < 0) ? previousEcg[previousEcg.Length + (i + del)] : ecgArray[i + del])
                                + a1 * ((i < 2 && hasEnoughDataFromPrevious && (i + del - 1) < 0) ? previousEcg[previousEcg.Length + (i + del - 1)] : ecgArray[i + del - 1])
                                + a2 * ((i < 2 && hasEnoughDataFromPrevious && (i + del - 2) < 0) ? previousEcg[previousEcg.Length + (i + del - 2)] : ecgArray[i + del - 2])
                                + b1 * ((hasEnoughDataFromPrevious && (i - 1) < 0) ? previousEcgFAMA[previousEcgFAMA.Length + (i - 1)] : ecgF[i - 1])
                                + b2 * ((hasEnoughDataFromPrevious && (i - 2) < 0) ? previousEcgFAMA[previousEcgFAMA.Length + (i - 2)] : ecgF[i - 2]);
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Crash: " + e.Message);
                        ecgF[i] = ecgArray[i];
                    }
                }

                return ecgF;

            }
            catch (Exception e)
            {
                Console.WriteLine("Crash: " + e.Message);
                return ecgArray;
            }
        }
    }