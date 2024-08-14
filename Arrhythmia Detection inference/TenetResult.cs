public class TenetResult : IEnumerable<double>
{
    public TenetResult(int[] peaks, bool[] nonECG, double[] noiseMeasure)
    {
        this.PeakIndices = peaks;
        this.NonECGInfo = nonECG;
        this.NoiseMeasure = noiseMeasure;
    }

    public int[] PeakIndices { get; set; }
    public bool[] NonECGInfo { get; set; }
    public double[] NoiseMeasure { get; set; }

    public IEnumerator<double> GetEnumerator()
    {
        return ((IEnumerable<double>)NoiseMeasure).GetEnumerator();
    }

    System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
    {
        return ((IEnumerable<double>)NoiseMeasure).GetEnumerator();
    }
}