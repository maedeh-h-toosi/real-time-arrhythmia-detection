using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

    public class Utils
    {

        public static T AvgIgnoreFake<T>(T[] array) where T : IComparable<T>
        {
            float res = 0;
            int FAKE_ECG_SAMPLE_VALUE = 40;
            if (array == null || array.Length == 0)
                return (T)(object)res;

            T fakeT = (T)Convert.ChangeType(FAKE_ECG_SAMPLE_VALUE, typeof(T));
            int count = 0;
            for (int j = 0; j < array.Length; j++)
            {
                T val = array[j];
                if (fakeT.CompareTo(val)!=0)
                {
                    res += Convert.ToSingle(val);
                    count += 1;
                }
            }
            res /= count;
            return (T)(object)res;
        }

    }
