using TorchSharp;
using static TorchSharp.torch;
using static TorchSharp.TensorExtensionMethods;
using System;
using Newtonsoft.Json;
using System.ComponentModel.DataAnnotations;
using ICSharpCode.SharpZipLib.Zip.Compression.Streams;
using System.IO;
using System.Globalization;
using ICSharpCode.SharpZipLib.Lzw;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
namespace FilterTenet;


class Program{
    
    static void Main(string[] args){


        int fs = 200;
        string csvFilePath = "C:/Users/tousi/Desktop/Inference_AF/first/new_test.csv";
        int desiredLength = 8000;
        List<float> ecgSignal = ReadECGFromCSV(csvFilePath);
        ecgSignal = PreProcess.EnsureDesiredLength(ecgSignal, desiredLength);
        int windowSize = 2000;
        List<List<float>> result = PreProcess.SlidingWindow(ecgSignal, windowSize);


        AdaptiveMedianFilter<float> filter = new AdaptiveMedianFilter<float>();
        AdaptiveMovingAverage butter = new AdaptiveMovingAverage();

        int[] rpeaks = new int[0];

        List<float[]> filteredECGSignals = new List<float[]>();
        List<(int index, string label)> segmentLabels = new List<(int, string)>();
        List<float[]> Filter_list = new List<float[]>();


        foreach (List<float> segment in result)
        {

            List<float> dividedSegment = segment.Select(x => (x - 2048) / 373f).ToList();
            float[] ecgData = dividedSegment.ToArray();
            

            float[] filteredECG = filter.Apply(ecgData, fs, rpeaks);
            float[] butterECG = butter.Apply_filter(filteredECG, fs);
            Filter_list.Add(butterECG);

        }

        double sourceSamplingRate = 200;
        double targetSamplingRate = 100;
        var model = new MobileNetV2("mobilenet");
        model.load("C:/Users/tousi/Desktop/Inference/first/new_mobilenet_model_weights.dat");
        model.eval();

        // Predict labels
        for (int i = 0; i < Filter_list.Count; i++)
        {
            var signal = Filter_list[i];
            double[] signalArray = signal.Select(x => (double)x).ToArray();
            double[] resampledSignal = PreProcess.CubicHermiteInterpol(signalArray, sourceSamplingRate, targetSamplingRate);
            List<float> resampledSignalList = resampledSignal.Select(x => (float)x).ToList();
            int f_desiredLength = 1000;
            resampledSignalList = PreProcess.EnsureDesiredLength(resampledSignalList, f_desiredLength);
            List<float> normalizedSignal = PreProcess.MinMax(resampledSignalList);

            Tensor inputTensor = torch.tensor(normalizedSignal.ToArray());
            Tensor outputTensor = model.forward(inputTensor);
            Tensor predictionTensor = outputTensor.argmax(1);
            long prediction = predictionTensor.item<long>();

            string label;
            switch (prediction)
            {
                case 0:
                    label = "Normal";
                    break;
                case 1:
                    label = "AF";
                    break;
                default:
                    label = $"Unknown class: {prediction}";
                    break;
            }
            segmentLabels[i] = (segmentLabels[i].index, label);
        }

        foreach (var segmentLabel in segmentLabels)
        {
            Console.WriteLine($"Segment {segmentLabel.index}: Predicted label: {segmentLabel.label}");
        }

     }
    



    
    static List<float> ReadECGFromCSV(string filePath)
    {
        List<float> ecgSignal = new List<float>();

        using (var reader = new StreamReader(filePath))
        {
            while (!reader.EndOfStream)
            {
                var line = reader.ReadLine();
                var values = line.Split(',').Select(float.Parse);
                ecgSignal.AddRange(values);
            }
        }

        return ecgSignal;
    }


    static void WriteResultsToCsv(string filePath, List<List<float>> results)
    {
        using (StreamWriter writer = new StreamWriter(filePath))
        {
            foreach (List<float> row in results)
            {
                string line = string.Join(",", row.Select(x => x.ToString())); 
                writer.WriteLine(line);
            }
        }
    }

}    
