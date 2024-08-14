using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Net.Mime.MediaTypeNames;
namespace FilterTenet;

    class QRSTenet
    {
        //Z-INTEGRATION: step 1
        const bool ENABLE_LOG_FILE = false;
        //Z-INTEGRATION: step 2
        /*static void Main(string[] args)
        {
            //Console.WriteLine(args[0]);
            //Console.WriteLine(args[1]);

            ///////////////////// Input arguments ///////////////////// 

            int fs = Int32.Parse(args[1]); //250;// 200; // 

            string filePath = args[0];  // "e0119.ekg";// "577556.ekg"; // //"232.ekg";  // replace with the path to your file

            /////////////////////////////////////////////////////
            */

        //Z-INTEGRATION: step 3
        public static TenetResult Main(float[] ecgSamples, int fs, out bool reachToEndSuccessfully)
        {
            //Z-INTEGRATION: step 4
            reachToEndSuccessfully = false;

            ///////////////////// Parameters ///////////////////// 

            ///////////////////////////////////// setting for non-ecg detector (SLAP)  begin
            int avg_segments = 1;
            int min_Necg_resolutions = 3;
            int min_ecg_resolutions = 20;
            Double comp = 0.6;  // maximum slap for non-ecg. if more than this, it is ecg
            Double hysp = 0.001;
            Double hysn = -0.001;
            Double Paus_Level = 0.05;

            Double half_slope_t = 0.01; ///  0.01 seconds
            Double half_resolution_t = 0.15; ///0.15 seconds
            int seg_to_res = 8; ///

            Double med_step = 0.3;  ///if median step is more than this, a non-ecg segment occurs
            ///////////////////////////////////// setting for non-ecg detector (SLAP)  end
            ///
            /// 
            ///////////////////////////////////// setting for QRSF and TENETSLF   begin
            Double win_t = 0.2;
            Double Bound_res_t = 0.005; // time resolution of Bound
            Double k_res_t = 0.025; // time resolution of k

            int peak_lim_portion = 6; // "peak limit portion" means the first 1/6 of peaks to be limited
            Double Decay_fin_ratio = 4.5;  // Decay final ratio 
            Double Decay_fin_lim = 0.55;  // Decay final limit 


            int DNF = 1; // Denoise Flag
            int TRF = 2; // T - Remove Flag
            Double QRSmin = 0.07;// Minimum voltage of QRS in mv;
            Double QRSmax = 6.0;// Maximum voltage of QRS in mv;
            int MagPow = 2;

            int cBound_min = 5;  // Minimum cBound;
            int cBound_max = 15; // Maximum cBound;

            // CBmin = 5; //% Change
            // CBmax = 15;
            //  CKmin = 3;
            // CKmax = 11; % Change

            int ck_min = 3; // Minimum ck;
            int ck_max = 11; // Maximum ck;

            ///////////////////////////////////// setting for    end
            ///////////////////////// setting for adaptive moving average
            Double f0min = 4;
            Double f0max = 40;
            Double Af = 0.9;
            Double Ar = 0.5;
            Double NFP = 0.8;
            Double NRP = 0.1;
            int delmax = 6;
            /////////////////////////////////////////////////////

            /*List<Double> lines = new List<Double>();

            // Read the file line by line and add each line to the list
            using (StreamReader reader = new StreamReader(filePath))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    lines.Add((Double)Convert.ToDouble(line));
                }
            }


            // Convert the list to an array
            Double[] linesArray = lines.ToArray();*/

            //Z-INTEGRATION: step 5
            Double[] PeakCurrent = null;
            int[] Seg = new int[ecgSamples.Length];
            try
            {
                //Z-INTEGRATION: step 6
                double[] linesArray = Array.ConvertAll(ecgSamples, x => (double)x);

                Double mean = Mean(linesArray, 0, linesArray.Length - 1);

                // Reduce all numbers from mean

                Double[] ekg = new Double[linesArray.Length];

                for (int i = 0; i < linesArray.Length; i++)
                {
                    ekg[i] = linesArray[i] - mean;
                }

                int length_ekg = linesArray.Length - 1;

                Double[] ekgr = new Double[length_ekg + 1];

                int Samples_in_half_slope = (int)(fs * half_slope_t);
                int Samples_in_half_resolution = (int)(fs * half_resolution_t);
                int Samples_in_half_segment = seg_to_res * Samples_in_half_resolution;

                //Z-INTEGRATION: step 9
                //int[] Seg = new int[length_ekg + 1];

                Double min_Necg = (2 * Samples_in_half_resolution + 1) * min_Necg_resolutions;

                /////////////////////////////////////

                ekgr[0] = ekg[0];

                ///////////////////////////////////// tries to remove large steps
                Double sstep = 0;

                for (int sample_number = 1; sample_number <= length_ekg; sample_number++)
                {
                    Double step = ekg[sample_number] - ekg[sample_number - 1];

                    if (Math.Abs(step) > (Double)(650 / fs))
                    {
                        sstep = sstep + step;
                        for (int i = sample_number; i <= Math.Min(sample_number + min_Necg - 1, length_ekg); i++)
                        {
                            Seg[i] = 1;
                        }

                    }
                    ekgr[sample_number] = ekg[sample_number] - sstep;

                }

                ///////////////////////////////////// calculate median from -1 second to 1 second 

                Double[] ecg = new Double[length_ekg + 1];
                ecg = MedFiltLL(ekgr, fs);


                for (int sample_number = fs; sample_number <= length_ekg - fs; sample_number++)
                {
                    if (Math.Abs(ecg[sample_number] - ecg[sample_number - 1]) > med_step)
                    {
                        for (int i1 = sample_number; i1 <= Math.Min(sample_number + min_Necg - 1, length_ekg); i1++)
                        {
                            Seg[i1] = 1;
                        }
                    }
                }

                for (int i = 0; i <= length_ekg; i++)
                    ecg[i] = ekgr[i] - ecg[i];

                //File.WriteAllLines("ecg_C.txt", Array.ConvertAll(ecg, x => x.ToString()));
                //File.WriteAllLines("ekgr_C.txt", Array.ConvertAll(ekgr, x => x.ToString()));
                //File.WriteAllLines("sstep_array_C.txt", Array.ConvertAll(sstep_array, x => x.ToString()));
                ///  File.WriteAllLines("ekg_C.txt", Array.ConvertAll(ekg, x => x.ToString()));
                ///////////////////////////////////// 

                ////////////////////////////////////// Using a filter var_w_2nd_adaptive_deq
                ///

                //List<Double> ecg,
                ///int fs; 


                List<Double[]> filterOut = var_w_2nd_adaptive_deq(ecg, fs, f0min, f0max, Af, Ar, NFP, NRP, delmax);
				Double[] ecgF = filterOut[0];
				Double[] noiseMeasure = filterOut[1];


                /////////////////////////////////////


                ///////////////////////////////////// non - ecg detector function --SLAPSeg function
                //var return_value_SLAPmag = SLAPSeg(ecg, Seg, Samples_in_half_slope, Samples_in_half_segment, Samples_in_half_resolution, avg_segments, min_Necg_resolutions, min_ecg_resolutions, comp, hysp, hysn, Paus_Level);
                var return_value_SLAPmag = SLAPSeg(ecgF, Seg, Samples_in_half_slope, Samples_in_half_segment, Samples_in_half_resolution, avg_segments, min_Necg_resolutions, min_ecg_resolutions, comp, hysp, hysn, Paus_Level);

                int WindowSize = (int)(win_t * fs);


                int ck_value = 0;
                int cBound_value = 0;

                //Z-INTEGRATION: step 10
                //Double[] PeakCurrent = null;
                Double[] slope = new Double[length_ekg + 1];
                Double avg = 0;
                Double[] avgA = new Double[length_ekg + 1];
                Double[] SLAP = new Double[length_ekg + 1];
                int[] Necg = new int[length_ekg + 1];
                int[] Paus = new int[length_ekg + 1];


                /*    File.WriteAllLines("SLAP.txt", Array.ConvertAll(SLAP, x => x.ToString()));
                    File.WriteAllLines("Necg.txt", Array.ConvertAll(Necg, x => x.ToString()));
                    File.WriteAllLines("Paus.txt", Array.ConvertAll(Paus, x => x.ToString()));
                    File.WriteAllLines("slope.txt", Array.ConvertAll(slope, x => x.ToString()));
                */
                Double[] High = new Double[length_ekg];

                slope = return_value_SLAPmag.Item1;
                avg = return_value_SLAPmag.Item2;
                avgA = return_value_SLAPmag.Item3;
                SLAP = return_value_SLAPmag.Item4;
                Necg = return_value_SLAPmag.Item5;
                Paus = return_value_SLAPmag.Item6;

                //  File.WriteAllLines("Necg.txt", Array.ConvertAll(Necg, x => x.ToString()));

                for (int i = 0; i <= length_ekg; i++)
                {
                    Seg[i] = Necg[i];
                }

                // Now Denosing is not working so ecgF = ecg


                /////////////////////
                ///  added for QRSF2

                Double[] SL = shiftRight(ecgF);        // Change
                SL[0] = ecgF[0];                // Change
                Double[] SR = shiftLeft(ecgF);       // Change
                SR[length_ekg] = ecgF[length_ekg]; // Change
                Double[] DL = new Double[length_ekg + 1];// zeros(size(ecgF));            // Change
                Double[] DR = new Double[length_ekg + 1]; ;                            // Change

                for (int i = 0; i < DL.Length; i++)
                {
                    DL[i] = 0;
                    DR[i] = 0;
                }

                int BoundC = 1;                         // Change

                List<Double> pDC = new List<Double>();
                List<int> pIDC = new List<int>();

                Double[] DMM = new Double[length_ekg + 1];
                Double[] BM = new Double[length_ekg + 1];

                /////////////////////////////////

                Double vte_Last = 0;
                for (int cBound = cBound_min; cBound <= cBound_max; cBound++)
                {
                    int Bound = (int)(cBound * fs * Bound_res_t);

                    Double[] DM = QRSF2(ecgF, Bound, DNF, TRF, QRSmin, QRSmax, MagPow, ref SL, ref SR, ref DL, ref DR, ref BoundC);


                    ////////////////////
                    ///

                    // new PVC det:
                    if (cBound == cBound_min)
                    {
                        for (int i = 0; i <= length_ekg; i++)
                        {
                            DMM[i] = DM[i];
                            BM[i] = Bound;
                        }

                    }
                    else
                    {
                        for (int sn = 0; sn <= length_ekg; sn++)
                        {
                            if (DM[sn] > DMM[sn])
                            {
                                DMM[sn] = DM[sn];
                                BM[sn] = Bound;
                            }
                        }
                    }
                    for (int i = 0; i <= length_ekg; i++)
                    {
                        DM[i] = 0.8 * DM[i] + 0.2 * DMM[i];
                    }

                    // new PVC det.


                    ////////////////////////////



                    (int[] pIDM, Double[] pDM) = findQRSpeakseg(DM, WindowSize);

                    /////  extraction of larger magnifier peaks   (possible R - peaks)
                    int LpDM = pIDM.Length - 1;
                    List<int> pIDMSI = new List<int>();


                    if (LpDM >= 0)
                        pIDMSI.Add(0);

                    int i2;

                    for (i2 = 1; i2 <= LpDM; i2++)
                    {
                        if (pDM[i2] >= pDM[pIDMSI[0]])
                            pIDMSI.Insert(0, i2);
                        else if (pDM[i2] <= pDM[pIDMSI[pIDMSI.Count - 1]])
                            pIDMSI.Add(i2);
                        else
                        {
                            int b4 = 0;
                            int b5 = pIDMSI.Count - 1;
                            int inc = 0;
                            while (b4 < b5 - 1)
                            {
                                int b6 = (int)((b4 + b5) / 2);
                                if (pDM[i2] <= pDM[pIDMSI[b6]])
                                {
                                    b4 = b6;
                                    inc = 1;
                                }
                                else
                                {
                                    inc = 0;
                                }
                                if (pDM[i2] >= pDM[pIDMSI[b6 + 1]])
                                {
                                    b5 = b6 + inc;
                                }
                                else
                                {
                                    b4 = b6 + 1;
                                }
                            }
                            pIDMSI.Insert(b4 + 1, i2);
                        }
                    }

                    int j = (int)Math.Ceiling((double)pIDMSI.Count / peak_lim_portion) - 1;
                    for (int i3 = 0; i3 < j; i3++)
                    {
                        pDM[pIDMSI[i3]] = pDM[pIDMSI[j]];
                    }

                    ////////////////////////////////////
                    int opt = 0;
                    for (int ck = ck_min; ck <= ck_max; ck++)
                    {
                        int k = (int)(ck * fs * k_res_t);
                        int SegStart = 0;
                        int pIDMSt = 0;
                        int pIDMSp = -1;
                        Double vte = 0;
                        int N = 0;   // N is the number of 

                        //List<Double> pDC = new List<Double>();
                        //List<int> pIDC = new List<int>();

                        pDC.Clear();
                        pIDC.Clear();


                        for (int sample_number = 1; sample_number <= length_ekg; sample_number++)
                        {
                            if ((pIDMSp < LpDM))                                    // new magnifier peak in the current segment
                            {
                                if ((pIDM[pIDMSp + 1] == sample_number - 1) || (sample_number == length_ekg))
                                    pIDMSp = pIDMSp + 1;
                            }

                            if ((sample_number == length_ekg) || !(Seg[sample_number] == Seg[sample_number - 1])) /// the condition for segment stop 
                            {
                                int[] pIDMSeg = new int[pIDMSp - pIDMSt + 1];

                                if (pIDMSp >= pIDMSt)                                           // creat pIDMSeg from pIDM
                                {
                                    Array.Copy(pIDM, pIDMSt, pIDMSeg, 0, pIDMSp - pIDMSt + 1);   //pIDMSeg = pIDM(pIDMSt: pIDMSp);

                                    for (int m = 0; m < pIDMSeg.Length; m++)                    //pIDMSeg = pIDMSeg - (SegStart - 1) * ones(size(pIDMSeg));
                                    {
                                        pIDMSeg[m] -= SegStart;
                                    }
                                }

                                //////////// TENETseg function implements decay in forward and backward time
                                Double[] pDM_tenet = new Double[pIDMSp - pIDMSt + 1];
                                Array.Copy(pDM, pIDMSt, pDM_tenet, 0, pIDMSp - pIDMSt + 1);

                                (int[] pIDCSeg, Double[] pDCSeg) = TENETSLFC(pDM_tenet, pIDMSeg, k, k, Math.Pow(Decay_fin_ratio, MagPow), Math.Pow(Decay_fin_lim, MagPow), WindowSize); ;     // Tenet function      
                                                                                                                                                                                              ///// decay to Vf which is a fraction of last peak  power(4.5, MagPow) and
                                if (pDCSeg != null)
                                    pDC.InsertRange(pDC.Count, pDCSeg);


                                if (pIDCSeg != null)
                                {
                                    for (int m = 0; m < pIDCSeg.Length; m++)
                                    {
                                        pIDCSeg[m] += SegStart;
                                    }
                                    pIDC.InsertRange(pIDC.Count, pIDCSeg);
                                }


                                /////////////////// consider ecg segment in ecg optimization

                                if (sample_number > 0)                                          //// we are in a ecg segment ((Seg(sample_number - 1) == 0)
                                {
                                    if (Seg[sample_number - 1] == 0)
                                    {
                                        int[] RRintervals = Diff(pIDCSeg);                        /// calculate RRintervals of current segment
                                        if (RRintervals != null)
                                        {
                                            int LRR = RRintervals.Length;
                                            if (LRR > 1)                                          // if there is more than RRintervals
                                            {
                                                N = N + LRR;                                    // N : number of valid RRintervals
                                                vte = vte + (LRR - 1) * Var(RRintervals);      // (LRR - 1) because of variance definiation
                                            }
                                        }
                                    }
                                }
                                SegStart = sample_number;                               // the next segment starts from this sample number 
                                pIDMSt = pIDMSp + 1;
                            }

                        } /// end for sample_number

                        if (N > 1)
                        {
                            vte = vte / (N - 1);
                        }

                        if (vte_Last == 0 || vte < vte_Last)   // the goal of optimization is to see a regular pattern 
                        {
                            cBound_value = cBound;
                            ck_value = ck;
                            vte_Last = vte;

                            opt = 1;

                            PeakCurrent = new Double[pIDC.Count];
                            Array.Copy(pIDC.ToArray(), PeakCurrent, pIDC.Count);          // PeakCurrent = pIDC;

                        }

                    }/// end for ck

                    ///// for show 
                    if (opt == 1)
                    {
                        /*    string[] str = new string[3];
                            str[0] = ck_value.ToString();
                            str[1] = cBound_value.ToString();
                            str[2] = vte_Last.ToString();
                            File.WriteAllLines(filePath + ".txt", str);


                            Console.WriteLine("cBound_value:" + cBound_value + "    ck_value: " + ck_value + "    vte_Last: " + vte_Last);
                        */
                        opt = 0;
                    }


                }/// end for cBound


                // new PVC det:
                Double[] Rmax = new Double[PeakCurrent.Length]; ////???????????
                Double[] Nmax = new Double[PeakCurrent.Length]; ////???????????
                Double[] SNmax = new Double[PeakCurrent.Length]; ////???????????
                for (int i = 0; i < PeakCurrent.Length; i++)  ///???????????
                {
                    Rmax[i] = 0;
                    Nmax[i] = 0;
                    SNmax[i] = 0;
                }

                int TW = (int)(fs * win_t / 3);
                for (int pn = 0; pn < PeakCurrent.Length; pn++) ///???????????///???????????
                {
                    for (int sn = (int)Math.Max(0, PeakCurrent[pn] - TW); sn < Math.Min(length_ekg, PeakCurrent[pn] + TW); sn++) ///???????????///???????????
                    {
                        if (DMM[sn] > Rmax[pn])
                        {
                            Rmax[pn] = DMM[sn];
                            Nmax[pn] = BM[sn];
                            SNmax[pn] = sn;
                        }
                    }
                }


                if (ENABLE_LOG_FILE) File.WriteAllLines("PeakCurrent_C#_SLFC4.txt", Array.ConvertAll(SNmax, x => x.ToString()));


                /// for debug
                /*
                string[] str = new string[2];
                str[0] = ck_value.ToString();
                str[1] = cBound_value.ToString();
                File.WriteAllLines(filePath+".txt", str);
                */

                //Z-INTEGRATION: step 7
                TenetResult result = new TenetResult(PeakCurrent.Select(p => (int)p).ToArray(), Seg.Select(s => s == 1).ToArray(), noiseMeasure);
                reachToEndSuccessfully = true;
                //Console.WriteLine("PeakIndices: " + string.Join(", ", result.PeakIndices));

                // Print NonECGInfo
               // Console.WriteLine("NonECGInfo: " + string.Join(", ", result.NonECGInfo));

                return result;
            }
            catch (Exception e)
            {
                //Z-INTEGRATION: step 8
                reachToEndSuccessfully = false;
                //Tools.SaveCatchLog(e, "Tenet Bugs...");
                int[] rpeaks = PeakCurrent == null ? new int[0] : PeakCurrent.Select(p => (int)p).ToArray();
                bool[] nonECGs = Seg == null ? new bool[ecgSamples.Length] : Seg.Select(s => s == 1).ToArray();
                double[] noiseMeasure = new double[0];
                TenetResult result = new TenetResult(rpeaks, nonECGs, noiseMeasure);
                Console.WriteLine(result);
                return result;
            }
        }
        static Double Var(int[] arr)
        {
            Double mean = arr.Average();
            Double sumOfSquares = (Double)arr.Sum(val => Math.Pow(val - mean, 2));
            return (Double)sumOfSquares / (arr.Length - 1);
        }
        static int[] Diff(int[] arr)
        {
            if (arr == null || arr.Length < 2)
                return null;

            if (arr.Length == 0)
            {
                Console.Write("arr.Length == 0 ");
            }
            int[] result = new int[arr.Length - 1];
            for (int i = 0; i < arr.Length - 1; i++)
            {
                result[i] = arr[i + 1] - arr[i];
            }
            return result;
        }

        static Double Mean(Double[] numbers, int start, int end)
        {
            if (numbers == null || numbers.Length == 0)
            {
                throw new ArgumentException("The array must not be null or empty.");
            }


            //Double sum = numbers.Skip(start).Take(end).Sum(); // returns 5.0
            Double sum = 0;
            for (int i = start; i <= end; i++)
            {
                sum += numbers[i];
            }

            // Double sum = numbers.Sum();
            Double mean = (Double)sum / (end - start + 1); // numbers.Length;

            return mean;
        }


        /// function  [pIDC, pDC] = TENETSLFC(pDM, pIDM, kF, kB, Vfr, Vfmax, WindowSize)
        public static (int[], Double[]) TENETSLFC(Double[] pDM, int[] pIDM, int kF, int kB, Double Vfr, Double Vfmax, int WindowSize)
        {

            int LpIM = pIDM.Length - 1;
            // Double[] pDC = new Double[LpIM + 1];
            List<Double> pDC_list = new List<Double>();
            List<int> pIDC_list = new List<int>();
            int[] pIDC;
            Double[] pDC;
            int i;

            if (LpIM < 1)
            {
                /* if (LpIM == -1)
                     return (null, null);
                 else {
                    // pDC[0] = pDM[0];
                     pDC_list.Add(pDM[0]);
                     pDC = pDC_list.ToArray();

                     pIDC_list.Add(pIDM[0]);
                     pIDC = pIDC_list.ToArray();

                     return ( pIDC, pDC);
              
                    }*/
                return (pIDM, pDM);
            }
            else
            {
                List<int> pII = new List<int>();
                int LpII = -1;
                if (LpIM > 1)
                {
                    int pC = 1;
                    for (i = 1; i <= (LpIM - 1); i++)
                    {

                        int cond1 = pC / 2;
                        int cond2 = (pDM[i] >= pDM[i - 1]) ? 1 : 0;
                        int cond3 = (pDM[i] >= pDM[i + 1]) ? 1 : 0;
                        if (cond1 + cond2 + cond3 >= 2)
                        {
                            pII.Add(i);
                            LpII = LpII + 1;
                            pC = 0;
                        }
                        else
                        {
                            pC = pC + 1;
                        }
                    }
                }

                //pDC[0] =pDM[0];
                pDC_list.Add(pDM[0]);
                Double Vf = Math.Min(Vfmax, pDC_list[0] / Vfr);
                i = 0;

                for (int peak_number = 1; peak_number <= LpIM; peak_number++)
                {
                    pDC_list.Add(Vf + (pDC_list[peak_number - 1] - Vf) * Math.Pow((Double)(kF - 1) / kF, pIDM[peak_number] - pIDM[peak_number - 1]));

                    if (pDM[peak_number] >= pDC_list[peak_number])
                    {
                        pDC_list[peak_number] = pDM[peak_number];
                        Vf = Math.Min(Vfmax, pDC_list[peak_number] / Vfr);
                    }
                    if (i <= LpII)
                    {
                        if (pII[i] == peak_number)
                        {
                            Vf = Math.Min(Vfmax, pDC_list[peak_number] / Vfr);
                            i = i + 1;
                        }
                    }
                }

                Vf = Math.Min(Vfmax, pDC_list[LpIM] / Vfr);
                i = LpII;

                for (int Reverse_peak_number = 1; Reverse_peak_number <= LpIM; Reverse_peak_number++)
                {
                    int peak_number = LpIM - Reverse_peak_number;///
                    Double DCS = Vf + (pDC_list[peak_number + 1] - Vf) * Math.Pow((Double)(kB - 1) / kB, pIDM[peak_number + 1] - pIDM[peak_number]);

                    if (pDC_list[peak_number] >= DCS)
                        Vf = Math.Min(Vfmax, pDC_list[peak_number] / Vfr);
                    else
                        pDC_list[peak_number] = DCS;

                    if (i >= 0)
                    {
                        if (pII[i] == peak_number)
                        {
                            Vf = Math.Min(Vfmax, pDC_list[peak_number] / Vfr);
                            i = i - 1;
                        }
                    }

                }

                int CondP = -1; // Change
                for (int peak_number = 0; peak_number <= LpIM; peak_number++)
                {

                    int Cond = (pDC_list[peak_number] == pDM[peak_number]) ? 1 : 0;

                    if ((Cond == 1) && (peak_number > 0))
                    {

                        Double preCond = pDC_list[peak_number - 1] * Math.Pow((Double)(kF - 1) / kF, pIDM[peak_number] - pIDM[peak_number - 1] - WindowSize);
                        Cond = Cond * ((pDC_list[peak_number] >= preCond) ? 1 : 0);
                    }
                    if ((Cond == 1) && (peak_number < LpIM))
                    {
                        Double preCond = pDC_list[peak_number + 1] * Math.Pow((Double)(kB - 1) / kB, pIDM[peak_number + 1] - pIDM[peak_number] - WindowSize);
                        Cond = Cond * ((pDC_list[peak_number] >= preCond) ? 1 : 0);
                    }
                    // if (Cond == 1)
                    //     pIDC_list.Add(pIDM[peak_number]);


                    ////////////////////
                    ///
                    if (Cond == 1)
                    {
                        if ((CondP == 0) && (pDM[Math.Max(0, peak_number - 1)] >= 0.7 * pDM[peak_number]) && (pIDM[peak_number] - pIDM[Math.Max(0, peak_number - 2)] < pIDM[Math.Min(LpIM, peak_number + 1)] - pIDM[peak_number]))  // Change
                        {
                            // LpIC = LpIC + 1; // Change
                            // pIDC[LpIC] = pIDM[peak_number - 1]; // Change
                            pIDC_list.Add(pIDM[peak_number - 1]);
                        }
                        // LpIC = LpIC + 1;
                        // pIDC[LpIC] = pIDM[peak_number];
                        pIDC_list.Add(pIDM[peak_number]);
                    }
                    CondP = Cond; // Change

                    ////////////

                }

            }



            pIDC = pIDC_list.ToArray();
            pDC = pDC_list.ToArray();

            /*  if (pIDC.Length != pDC.Length)
              {
                  Console.WriteLine("  pIDC.Length != pDC.Length  in TENETSLFC function ");
              }
            */
            return (pIDC, pDC);
        }

        public static (Double, int) FindMax(Double[] arr, int start, int stop)
        {
            int maxIndex = start;
            Double maxValue = arr[start];
            for (int i = start + 1; i <= stop; i++)
            {
                if (arr[i] > maxValue)
                {
                    maxIndex = i;
                    maxValue = arr[i];
                }
            }
            maxIndex = maxIndex - start;
            return (maxValue, maxIndex);
        }


        public static (Double, int) FindMin(Double[] arr, int start, int stop)
        {
            int minIndex = start;
            Double minValue = arr[start];
            for (int i = start + 1; i <= stop; i++)
            {
                if (arr[i] < minValue)
                {
                    minIndex = i;
                    minValue = arr[i];
                }
            }
            minIndex = minIndex - start;
            return (minValue, minIndex);
        }
        //function  [slope, avg, avgA, SLAP, Necg, Paus] = SLAPmag(ecg, Samples_in_half_slope, Samples_in_half_segment, Samples_in_half_resolution, avg_segments, min_Necg_resolutions, min_ecg_resolutions, comp, hysp, hysn, Paus_Level)

        public static Tuple<Double[], Double, Double[], Double[], int[], int[]> SLAPSeg(Double[] ecg, int[] Seg, int Samples_in_half_slope,
                                                                            int Samples_in_half_segment, int Samples_in_half_resolution,
                                                                            int avg_segments, int min_Necg_resolutions, int min_ecg_resolutions,
                                                                            Double comp, Double hysp, Double hysn, Double Paus_Level)
        {

            int length_ekg = ecg.Length - 1;


            Double[] slope = new Double[length_ekg + 1];
            Double avg = 0;
            Double[] avgA = new Double[length_ekg + 1];
            Double[] SLAP = new Double[length_ekg + 1];
            int[] Necg = new int[length_ekg + 1];
            int[] Paus = new int[length_ekg + 1];

            Double[] High = new Double[length_ekg + 1];
            Double hyst = hysn;

            for (int i = 0; i <= length_ekg; i++)
            {
                slope[i] = 0;
                SLAP[i] = 0;
                High[i] = 0;
                Necg[i] = 0;
                Paus[i] = 0;
                avgA[i] = 0;
            }

            int st;
            int sp;

            for (int sample_number = 0; sample_number <= length_ekg; sample_number++)
            {
                st = Math.Max(sample_number - Samples_in_half_slope, 0);
                sp = Math.Min(sample_number + Samples_in_half_slope, length_ekg);
                slope[sample_number] = 0;  //
                for (int i = st; i <= sp; i++)
                {
                    slope[sample_number] = slope[sample_number] + (2 * i - (sp + st)) * ecg[i];//
                }

                slope[sample_number] = Math.Abs(6 * slope[sample_number] / ((sp - st + 1) * (sp - st + 2)));//
            }


            avg = Mean(slope, 0, length_ekg);


            Double avgAi = 0;
            int avg_sp = 0;
            int avg_st = 0;
            int dsn_upperBound = (int)((length_ekg - 2 * Samples_in_half_segment) / (2 * Samples_in_half_resolution + 1));




            for (int dsn = 0; dsn <= dsn_upperBound; dsn++)
            {
                st = dsn * (2 * Samples_in_half_resolution + 1);
                sp = st + 2 * Samples_in_half_segment;
                int avg_st_n = Math.Max(0, Math.Min(st - (avg_segments - 1) * Samples_in_half_segment, length_ekg - avg_segments * (2 * Samples_in_half_segment + 1)));
                int avg_sp_n = avg_st_n + avg_segments * (2 * Samples_in_half_segment + 1);

                if (dsn == 0)
                    avgAi = Mean(slope, avg_st_n, avg_sp_n);
                else if ((avg_sp < avg_sp_n) && (avg_st < avg_st_n))
                {
                    // avgAi = avgAi + slope.Skip(avg_sp + 1).Take(avg_sp_n).Sum() - slope.Skip(avg_st).Take(avg_st_n - 1).Sum() / (avg_sp_n - avg_st_n + 1);
                    Double sum_sp = 0;
                    Double sum_st = 0;
                    for (int i = avg_sp + 1; i <= avg_sp_n; i++)
                    {
                        sum_sp += slope[i];
                    }
                    for (int i = avg_st; i <= avg_st_n - 1; i++)
                    {
                        sum_st += slope[i];
                    }
                    avgAi = avgAi + (sum_sp - sum_st) / (avg_sp_n - avg_st_n + 1);
                }
                avgA[st + Samples_in_half_segment] = avgAi;
                avg_st = avg_st_n;
                avg_sp = avg_sp_n;

                //  if( dsn % 100 == 0)
                //      Console.WriteLine(st + "   "+"   " + sp + "   " + avg_st_n + "   " + avg_sp_n + "   " + tmp);

                Double ekgmax = ecg[st];
                Double ekgmin = ekgmax;
                for (int i = st; i <= sp; i++)
                {
                    if (slope[i] < avgAi)
                    {
                        SLAP[st + Samples_in_half_segment] = SLAP[st + Samples_in_half_segment] + 1; // the number of slopes less than average (SLAP)  P refers to portion
                    }
                    if (ecg[i] > ekgmax)
                    {
                        ekgmax = ecg[i];
                    }
                    if (ecg[i] < ekgmin)
                    {
                        ekgmin = ecg[i];
                    }

                }


                SLAP[st + Samples_in_half_segment] = SLAP[st + Samples_in_half_segment] / (sp - st + 1);

                if ((ekgmax - ekgmin) > 10)  /// 10 means 10 mv
                {
                    for (int i = st; i <= sp; i++)
                    {
                        if (High[i] == 0)
                        {
                            High[i] = 1;
                        }
                    }
                }
                else
                {
                    for (int i = st; i <= sp; i++)
                    {
                        High[i] = -1;
                    }
                }

                if (SLAP[st + Samples_in_half_segment] < (comp + hyst))
                {
                    for (int i = st; i <= sp; i++)
                    {
                        Necg[i] = Necg[i] + 1;
                    }
                    hyst = hysp;
                }
                else
                {
                    for (int i = st; i <= sp; i++)
                    {
                        Necg[i] = Necg[i] - 1;
                    }
                    hyst = hysn;
                }

                if ((ekgmax - ekgmin) < Paus_Level)
                {
                    for (int i = st; i <= sp; i++)
                    {
                        Paus[i] = 1;
                    }
                }


            }


            for (int sample_number = 0; sample_number <= length_ekg; sample_number++)
            {
                if (Necg[sample_number] <= 0 && High[sample_number] <= 0 && Seg[sample_number] == 0)
                    Necg[sample_number] = 0;
                else
                    Necg[sample_number] = 1;
            }

            int stf = 0; // start flag
            Necg[length_ekg] = 1;
            sp = 0;
            for (int sample_number = 0; sample_number <= length_ekg; sample_number++)// remove short non-ecg segments 
            {
                if (Necg[sample_number] == 0 && stf == 0)
                {
                    stf = 1;
                    st = sample_number;
                    if ((st - sp) < (min_Necg_resolutions * (2 * Samples_in_half_resolution + 1)))
                    {
                        for (int sn = sp; sn <= st; sn++) //previous non-ecg segment
                        {
                            if (High[sn] <= 0)
                            {
                                Necg[sn] = 0;
                            }
                        }
                    }
                }

                if (Necg[sample_number] == 1 && stf == 1) // we arrive the end of current ecg segment
                {
                    stf = 0;
                    sp = sample_number;
                }
            }


            stf = 0;
            sp = 0;
            st = 1; ////???????????????????
            for (int sample_number = 0; sample_number <= length_ekg; sample_number++) // remove short ecg segments 
            {
                if (Necg[sample_number] == 0 && stf == 0)
                {
                    stf = 1;
                    st = sample_number;
                }
                if (Necg[sample_number] == 1 && stf == 1)
                {
                    stf = 0;
                    sp = sample_number;
                    if ((sp - st) < (min_ecg_resolutions * (2 * Samples_in_half_resolution + 1)))
                    {
                        for (int sn = st; sn <= sp; sn++)
                        {
                            Necg[sn] = 1;
                        }
                    }
                }

            }


            return Tuple.Create(slope, avg, avgA, SLAP, Necg, Paus);

        }
        public static Double[] QRSeg(Double[] ecg, int Bound, int DenoiseFlag, int TRemoveFlag, Double QRSmin, Double QRSmax, int MagPow)
        {
            int length = ecg.Length - 1;
            Double[] DM = new Double[length + 1];

            Double e_i;
            Double ei;
            Double e_i_n;
            Double ein;

            for (int sample_number = 0; sample_number <= length; sample_number++)
            {
                Double M = 0;
                Double DMi = 0;
                Double DM_i = 0;


                for (int i = 0; i <= Bound; i++)
                {
                    if (sample_number >= i)
                        e_i = ecg[sample_number - i];
                    else
                        e_i = 0;

                    if ((sample_number + i) <= length)
                        ei = ecg[sample_number + i];
                    else
                        ei = 0;

                    if (i == 0 && sample_number > 0 && sample_number < length)
                    {
                        ei = ecg[sample_number + 1];
                        e_i = ecg[sample_number - 1];
                    }

                    if (sample_number >= i + Bound)
                        e_i_n = ecg[sample_number - i - Bound];
                    else
                        e_i_n = 0;

                    if (sample_number + i + Bound <= length)
                        ein = ecg[sample_number + i + Bound];
                    else
                        ein = 0;


                    DMi = DMi + (Bound - 2 * i) * ei;
                    DM_i = DM_i + (Bound - 2 * i) * e_i;
                    DM[sample_number] = DMi + DM_i;
                    M = M + (Bound - 2 * i) * (e_i_n + ein);

                }

                if (M * DM[sample_number] < 0) // Noise condition
                {
                    if (DenoiseFlag == 1)
                        DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                    else if (DenoiseFlag != 0)
                    {
                        Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                        Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                        DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                    }
                }
                else                            // T-wave condition
                    if (TRemoveFlag == 1)
                    DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                else if (TRemoveFlag != 0)
                {
                    Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                    Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                    DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                }

                DM[sample_number] = 6 * (Double)Math.Abs(DM[sample_number]) / (Double)((Bound + 1) * (Bound + 2));
                DM[sample_number] = (Double)Math.Pow(Math.Min(Math.Max(DM[sample_number], 2 * QRSmin), 2 * QRSmax), MagPow);

                // DM[sample_number] = System.Math.Round(DM[sample_number], 7);   
            }

            return DM;
        }

        public static Double[] QRSF2(Double[] ecg, int Bound, int DenoiseFlag, int TRemoveFlag, Double QRSmin, Double QRSmax, int MagPow, ref Double[] SL, ref Double[] SR, ref Double[] DL, ref Double[] DR, ref int BoundC)
        {
            int Lecg = ecg.Length - 1;

            Double[] DM = new Double[Lecg + 1];  // Masoud changed length to length + 1

            //Console.WriteLine("Begin Bound: " + Bound);
            // Console.WriteLine("Begin BoundC: " + BoundC);

            for (BoundC = BoundC + 1; BoundC <= Bound; BoundC++)
            {
                for (int sn = 0; sn <= Lecg; sn++)
                {

                    SL[sn] = SL[sn] + ecg[Math.Max(0, sn - BoundC + 1)];
                    SR[sn] = SR[sn] + ecg[Math.Min(Lecg, sn + BoundC - 1)];
                    DL[sn] = DL[sn] + SL[sn] - BoundC * ecg[Math.Max(0, sn - BoundC)];
                    DR[sn] = DR[sn] + SR[sn] - BoundC * ecg[Math.Min(Lecg, sn + BoundC)];
                }
            }
            //Console.WriteLine("End Bound: " + Bound);
            // Console.WriteLine("End BoundC: " + BoundC);
            // System.Environment.Exit(0);
            BoundC = BoundC - 1;  // very importaant because the difference of matlab and c in for loop

            for (int i = 0; i <= Lecg; i++)
            {
                DM[i] = DL[i] + DR[i];
            }


            for (int sn = 0; sn <= Lecg; sn++)
            {
                Double M = DL[Math.Max(0, sn - Bound)] + DR[Math.Min(Lecg, sn + Bound)];

                if (M * DM[sn] < 0)
                {
                    if (DenoiseFlag == 1)
                        DM[sn] = DM[sn] - Math.Sign(DM[sn]) * Math.Min(Math.Abs(M), Math.Abs(DM[sn]));
                    else if (DenoiseFlag != 0)
                        DM[sn] = DM[sn] - Math.Sign(DM[sn]) * Math.Min(Math.Min(Math.Abs(M), Math.Abs(DL[sn])), Math.Min(Math.Abs(DR[sn]), Math.Abs(DM[sn])));
                }
                else
                {

                    if (TRemoveFlag == 1)
                        DM[sn] = DM[sn] - Math.Sign(DM[sn]) * Math.Min(Math.Abs(M), Math.Abs(DM[sn]));
                    else if (TRemoveFlag != 0)
                        DM[sn] = DM[sn] - Math.Sign(DM[sn]) * Math.Min(Math.Min(Math.Abs(M), Math.Abs(DL[sn])), Math.Min(Math.Abs(DR[sn]), Math.Abs(DM[sn])));
                }

                DM[sn] = 6 * Math.Abs(DM[sn]) / (Double)((Bound + 1) * (Bound + 2));

                DM[sn] = Math.Pow(Math.Min(Math.Max(DM[sn], 2 * QRSmin), 2 * QRSmax), MagPow);

            }

            return DM;
        }


        /*    public static Double[] QRSF(Double[] ecg, int Bound, int DenoiseFlag, int TRemoveFlag, Double QRSmin, Double QRSmax, int MagPow)
            {
                int length = ecg.Length - 1;
                Double[] DM = new Double[length + 1];  // Masoud changed length to length + 1

                Double e_i;
                Double ei;
                Double e_i_n;
                Double ein;

                for (int sample_number = 0; sample_number <= Math.Min(2 * Bound - 1, length); sample_number++)
                {
                    Double M = 0;
                    //DM[sample_number - 1] = 0;

                    Double DMi = 0;
                    Double DM_i = 0;

                    for (int i = 0; i <= Bound; i++)
                    {
                        if (sample_number >= i)
                            e_i = ecg[sample_number - i];  // sample_number - i   
                        else
                            e_i = 0;

                        if ((sample_number + i) <= length)   ///    
                            ei = ecg[sample_number + i];    ///    sample_number + i  
                        else
                            ei = 0;

                        if (i == 0 && sample_number > 0 && sample_number < length)
                        {
                            ei = ecg[sample_number + 1];     /// sample_number + 1      
                            e_i = ecg[sample_number - 1];   ///  sample_number - 1  
                        }

                        if (sample_number >= i + Bound)
                            e_i_n = ecg[sample_number - i - Bound]; //   sample_number - i - Bound    
                        else
                            e_i_n = 0;

                        if (sample_number + i + Bound <= length)   ///    
                            ein = ecg[sample_number + i + Bound];   /// sample_number + i + Bound
                        else
                            ein = 0;

                        // DM[sample_number - 1] = DM[sample_number - 1] + (Bound - 2 * i) * (e_i + ei);
                        //  M = M + (Bound - 2 * i) * (e_i_n + ein);

                        DMi = DMi + (Bound - 2 * i) * ei;
                        DM_i = DM_i + (Bound - 2 * i) * e_i;

                        M = M + (Bound - 2 * i) * (e_i_n + ein);

                    }
                    DM[sample_number] = DMi + DM_i;

                    if (M * DM[sample_number] < 0) // Noise condition
                    {
                        if (DenoiseFlag == 1)
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                        else if (DenoiseFlag != 0)
                        {
                            Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                            Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                        }
                    }
                    else
                    {                        // T-wave condition
                        if (TRemoveFlag == 1)
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                        else if (TRemoveFlag != 0)
                        {
                            Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                            Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                        }
                    }

                    DM[sample_number] = 6 * (Double)Math.Abs(DM[sample_number]) / (Double)((Bound + 1) * (Bound + 2));
                    DM[sample_number] = (Double)Math.Pow(Math.Min(Math.Max(DM[sample_number], 2 * QRSmin), 2 * QRSmax), MagPow);

                    // DM[sample_number] = System.Math.Round(DM[sample_number], 7);   
                }


                for (int sample_number = 2 * Bound; sample_number <= length - 2 * Bound; sample_number++)
                {
                    Double M = 0;
                    //DM[sample_number - 1] = 0;

                    Double DMi = 0;
                    Double DM_i = 0;

                    //   if (sample_number == 140260) { 
                    //       Console.WriteLine("test");
                    //    }
                    int C = Bound;

                    for (int i = 0; i <= Bound; i++)
                    {
                        DMi = DMi + C * ecg[sample_number + Math.Max(i, 1)];
                        DM_i = DM_i + C * ecg[sample_number - Math.Max(i, 1)];
                        M = M + C * (ecg[sample_number - i - Bound] + ecg[sample_number + i + Bound]);
                        C -= 2;

                    }
                    DM[sample_number] = DMi + DM_i;


                    if (M * DM[sample_number] < 0) // Noise condition
                    {
                        if (DenoiseFlag == 1)
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                        else if (DenoiseFlag != 0)
                        {
                            Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                            Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                        }
                    }
                    else
                    {                            // T-wave condition
                        if (TRemoveFlag == 1)
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                        else if (TRemoveFlag != 0)
                        {
                            Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                            Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                        }
                    }

                    DM[sample_number] = 6 * (Double)Math.Abs(DM[sample_number]) / (Double)((Bound + 1) * (Bound + 2));
                    DM[sample_number] = (Double)Math.Pow(Math.Min(Math.Max(DM[sample_number], 2 * QRSmin), 2 * QRSmax), MagPow);

                    // DM[sample_number] = System.Math.Round(DM[sample_number], 7);   
                }

                for (int sample_number = Math.Max(length - 2 * Bound + 1, 2 * Bound); sample_number <= length; sample_number++)
                {
                    Double M = 0;
                    Double DMi = 0;
                    Double DM_i = 0;


                    for (int i = 0; i <= Bound; i++)
                    {

                        e_i = ecg[sample_number - i];  // sample_number - i   

                        if ((sample_number + i) <= length)   ///    
                            ei = ecg[sample_number + i];    ///    sample_number + i  
                        else
                            ei = 0;

                        if (i == 0 && sample_number < length)
                        {
                            ei = ecg[sample_number + 1];     /// sample_number + 1      
                            e_i = ecg[sample_number - 1];   ///  sample_number - 1  
                        }

                        if (sample_number + i + Bound <= length)   ///    
                            ein = ecg[sample_number + i + Bound];   /// sample_number + i + Bound
                        else
                            ein = 0;


                        DMi = DMi + (Bound - 2 * i) * ei;
                        DM_i = DM_i + (Bound - 2 * i) * e_i;

                        M = M + (Bound - 2 * i) * (ecg[sample_number - i - Bound] + ein);

                    }
                    DM[sample_number] = DMi + DM_i;

                    if (M * DM[sample_number] < 0) // Noise condition
                    {
                        if (DenoiseFlag == 1)
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                        else if (DenoiseFlag != 0)
                        {
                            Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                            Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                        }
                    }
                    else
                    {                          // T-wave condition
                        if (TRemoveFlag == 1)
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(Math.Abs(M), Math.Abs(DM[sample_number]));
                        else if (TRemoveFlag != 0)
                        {
                            Double min1 = Math.Min(Math.Abs(M), Math.Abs(DMi));
                            Double min2 = Math.Min(Math.Abs(DM_i), Math.Abs(DM[sample_number]));
                            DM[sample_number] = DM[sample_number] - Math.Sign(DM[sample_number]) * Math.Min(min1, min2);
                        }
                    }
                    DM[sample_number] = 6 * (Double)Math.Abs(DM[sample_number]) / (Double)((Bound + 1) * (Bound + 2));
                    DM[sample_number] = (Double)Math.Pow(Math.Min(Math.Max(DM[sample_number], 2 * QRSmin), 2 * QRSmax), MagPow);

                    // DM[sample_number] = System.Math.Round(DM[sample_number], 7);   
                }

                return DM;
            }
        */
        public static (int[], Double[]) findQRSpeakseg(Double[] ecg, int win)
        {
            //int length = ;
            //Double[] peakI = new Double[length];
            //Double[] peakAmp = new Double[length];
            //  Console.WriteLine("length ecg is: " + length + " win " + win);

            List<int> list_peakI = new List<int>();
            List<Double> list_peakAmp = new List<Double>();

            int L = ecg.Length - 1;
            int start = 1;  // new modif
            //int ind = 0;

            while (start < L)
            {
                int istart = start;
                int stop = Math.Min(start + win - 1, L);
                int stopn = stop;
                Double M;
                int I;

                (M, I) = FindMax(ecg, start, stop);

                if (stop < L && I > 0) // new modif
                {
                    int startn = stop + 1;
                    Double M1;
                    int I1;
                    stopn = Math.Min(start + win + I - 1, L);  // new modif

                    (M1, I1) = FindMax(ecg, startn, stopn);

                    while (M1 > M)
                    {
                        M = M1;
                        I = I1;
                        start = startn;
                        stop = stopn;
                        if (stop < L)
                        {
                            startn = stop + 1;
                            stopn = Math.Min(start + win + I - 1, L);  // new modif

                            (M1, I1) = FindMax(ecg, startn, stopn);

                        }
                    }
                }

                int peI = start + I;  //int peI = start + I - 1;
                //Console.WriteLine(Math.Min(Math.Min(peI - win, 1), istart));

                // new modif  start 
                int cond = 0;

                // reject the last sampls as a peak
                if (peI == L)
                {
                    cond = 0;
                }
                else  // it is not the last sample and can be a peak
                {
                    cond = peI > istart ? 1 : 0;
                    for (int sn = Math.Max(peI - win + 1, 0); sn <= istart - 1; sn++)
                    {
                        if (ecg[sn] > M)
                        {
                            cond = 0;
                            break;
                        }
                        if (ecg[sn] < M)
                        {
                            cond = 1;
                        }
                    }

                }

                if (cond == 1)   // new modif  end
                {
                    list_peakI.Add(peI);
                    list_peakAmp.Add(M);
                }

                start = stopn + 1;


            }

            int[] Array_peakI = list_peakI.ToArray();
            Double[] Array_peakAmp = list_peakAmp.ToArray();
            return (Array_peakI, Array_peakAmp);
        }

        public static Double Sum(Double[] numbers, int start, int end)
        {
            Double sum = 0;
            for (int number = start; number < end; number++)
            {
                sum += number;
            }
            return sum;
        }

        public static int SearchSortedLL(List<Double> SortedArray, int ArrLen, Double Value) // Arrlen is the index of the last SortedArray elemets
        {
            int Sample_Num;
            if (Value >= SortedArray[0])
                Sample_Num = -1;
            else if (Value <= SortedArray[ArrLen])
                Sample_Num = ArrLen;
            else
            {
                int B4 = 0;
                int B5 = ArrLen;
                while (B4 < B5 - 1)
                {
                    int B6 = (int)((B4 + B5) / 2);
                    int B4Flag = 0;

                    if (Value <= SortedArray[B6])
                    {
                        B4 = B6;
                        B4Flag = 0;
                    }
                    else
                        B4Flag = 1;

                    if (Value >= SortedArray[B6 + 1])
                        B5 = B6 + 1 - B4Flag;
                    else
                        B4 = B6 + 1;

                }
                if (Value == SortedArray[B5])
                    Sample_Num = B5;
                else
                    Sample_Num = B4;
            }

            return Sample_Num;

        }

        public static Double[] MedFiltLL(Double[] Arr, int Half_Window)
        {
            int ArrLen = Arr.Length - 1;
            int WindowSize = 2 * Half_Window + 1;

            //Double[] SortedArray = new Double[WindowSize];// zeros(2 * Half_Window + 1, 1);
            Double[] Med = new Double[ArrLen + 1];  //zeros(ArrLen, 1);

            List<Double> SortedArray = new List<Double>();

            SortedArray.Add(Arr[0]); //SortedArray[0] = Arr[0];
            int Insert_Sorted_Num;
            int Remove_Sorted_Num;

            for (int Sample_Num = 1; Sample_Num <= ArrLen; Sample_Num++)
            {
                if (Sample_Num >= WindowSize)
                {
                    Insert_Sorted_Num = SearchSortedLL(SortedArray, WindowSize - 1, Arr[Sample_Num]);
                    Remove_Sorted_Num = Math.Max(0, SearchSortedLL(SortedArray, WindowSize - 1, Arr[Sample_Num - WindowSize]));
                }
                else
                {
                    Insert_Sorted_Num = SearchSortedLL(SortedArray, Sample_Num - 1, Arr[Sample_Num]);
                    Remove_Sorted_Num = Sample_Num;
                    SortedArray.Add(0);
                }

                /* int SN = Remove_Sorted_Num;
                 int Dir = Math.Sign(Insert_Sorted_Num - Remove_Sorted_Num);
                 int Fin = Insert_Sorted_Num - Math.Min(Dir, 0);
                 while (!(SN == Fin))
                 {
                     SortedArray[SN] = SortedArray[SN + Dir];
                     SN = SN + Dir;
                 }
                 SortedArray[SN] = Arr[Sample_Num];
                */
                if (Insert_Sorted_Num < Remove_Sorted_Num)
                {
                    SortedArray.RemoveAt(Remove_Sorted_Num);
                    SortedArray.Insert(Insert_Sorted_Num + 1, Arr[Sample_Num]);
                }
                else
                {
                    SortedArray.Insert(Insert_Sorted_Num + 1, Arr[Sample_Num]);
                    SortedArray.RemoveAt(Remove_Sorted_Num);
                }

                if (Sample_Num >= 2 * Half_Window)
                    Med[Sample_Num - Half_Window] = SortedArray[Half_Window];

            }

            for (int Sample_Num = 0; Sample_Num < Half_Window; Sample_Num++)
            {
                Med[Sample_Num] = Med[Half_Window];
                Med[ArrLen - Sample_Num] = Med[ArrLen - Half_Window];
            }

            return Med;
        }

        /*    public static Double MyPow(Double input , int MP)
            {
                Double output = 1;

                for (int i = 0; i < MP; i++)
                {
                    output = output * input;

                }

                return output;
            }
        */

        public static List<Double[]> var_w_2nd_adaptive_deq(Double[] ecg, int fs, Double f0min, Double f0max, Double Af, Double Ar, Double NFP, Double NRP, int delmax)
        {
            //  fs:  sampling frequency
            //  Af: higher (0 < A_fall < 1) for slower cut-off fall
            //  Ar: lower (0 < A_rise < 1) for faster cut-off rise
            //  f0min  min cut - off frequency
            //  f0max  cut-off frequency


            Double p = Math.Sqrt(2);  //Critically Damped : 2, Bessel : 3, Butterworth : sqrt(2) 
            Double g = 1;  //Critically Damped : 1, Bessel : 3, Butterworth : 1
            Double c = Math.Sqrt(2 / (2 * g - p * p + Math.Sqrt(Math.Pow(2 * g - p * p, 2) + 4 * g * g)));


            int Lecg = ecg.Length;

            Double[] ecgF = new Double[Lecg];
            for (int i = 0; i < Lecg; i++)
            {
                ecgF[i] = ecg[i];
            }

            Double f0 = f0max;
            Double NM = 1.5;
			List<Double> NoiseMeasure = new List<Double> { };

            for (int i = 2; i < Lecg; i++)
            {
                Double fc = 0;

                for (int j = 0; j <= delmax; j++)
                {

                    fc = fc + 120 * f0min * Math.Abs(ecg[Math.Max(0, Math.Min(i + j, Lecg - 1))] - ecgF[i - 1]) / (Double)(delmax + 1);
                }

                fc = Math.Max(Math.Max(f0min, f0 * Af), fc / NM);
                f0 = Math.Min(Math.Min(f0max, f0 / Ar), fc);

                NM = Math.Max(0.3, NM * (1 + Math.Min(NRP, (f0 / f0min) - 1 - NFP) / 100));
				NoiseMeasure.Add(NM);

                Double pi = 3.1415926;
                Double w0 = Math.Tan(pi * c * f0 / fs); //digital cut - off frequency
                Double a0 = g * w0 * w0 / (1 + p * w0 + g * w0 * w0);
                Double a1 = 2 * a0;
                Double a2 = a0;
                Double b1 = 2 * a0 * (1 / (g * w0 * w0) - 1);
                Double b2 = 1 - a0 - a1 - a2 - b1;

                // ecgF[i] = a0 * ecg[i] + a1 * ecg[i - 1] + a2 * ecg[i - 2] + b1 * ecgF[i - 1] + b2 * ecgF[i - 2];

                int del = (int)Math.Min(Math.Min(delmax, Math.Round(0.24 * fs / f0)), Lecg - 1 - i);
                ecgF[i] = a0 * ecg[i + del] + a1 * ecg[i + del - 1] + a2 * ecg[i + del - 2] + b1 * ecgF[i - 1] + b2 * ecgF[i - 2];

            }

           

            //File.WriteAllLines("ecgF_var_moving_avg_C_EQualized.txt", Array.ConvertAll(ecgF, x => x.ToString()));
            return new List<Double[]> { ecgF, NoiseMeasure.ToArray() };
        }

        public static Double[] shiftLeft(Double[] arr)
        {
            var result = new Double[arr.Length];
            Array.Copy(arr, 1, result, 0, arr.Length - 1);
            result[arr.Length - 1] = arr[0];
            return result;

        }

        public static Double[] shiftRight(Double[] arr)
        {
            var result = new Double[arr.Length];
            Array.Copy(arr, 0, result, 1, arr.Length - 1);
            result[0] = arr[arr.Length - 1];
            return result;
        }
    }