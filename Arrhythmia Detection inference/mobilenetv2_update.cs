// Copyright (c) .NET Foundation and Contributors.  All Rights Reserved.  See LICENSE in the project root for license information.

// A number of implementation details in this file have been translated from the Python version of torchvision,
// largely located in the files found in this folder:
//
// https://github.com/pytorch/vision/blob/a4f53308b2d0f1aa9191686e326f45c26053f686/torchvision/models/mobilenetv2.py
//
// The origin has the following copyright notice and license:
//
// https://github.com/pytorch/vision/blob/main/LICENSE
//

using System;
using System.Collections.Generic;

using static TorchSharp.torch;

// namespace TorchSharp
// {
//     namespace Modules
//     {
        /// <summary>
        /// MobileNet V2 main class
        /// </summary>
        public class MobileNetV2 : nn.Module<Tensor, Tensor>
        {
        
            private class InvertedResidual : nn.Module<Tensor, Tensor>
            {
                private readonly bool _is_cn;
                private readonly nn.Module<Tensor, Tensor> conv;
                private readonly long out_channels;
                private readonly long stride;
                private readonly bool use_res_connect;

                protected override void Dispose(bool disposing)
                {
                    if (disposing) {
                        conv.Dispose();
                    }
                    base.Dispose(disposing);
                }

                public InvertedResidual(
                    string name,
                    long inp,
                    long oup,
                    long stride,
                    double expand_ratio,
                    Func<long, nn.Module<Tensor, Tensor>>? norm_layer = null) : base(name)
                {
                    this.stride = stride;
                    if (stride != 1 && stride != 2) {
                        throw new ArgumentOutOfRangeException($"stride should be 1 or 2 insted of {stride}");
                    }

                    if (norm_layer == null) {
                        norm_layer = (features) => nn.BatchNorm1d(features);
                    }

                    var hidden_dim = (long)Math.Round(inp * expand_ratio);
                    this.use_res_connect = this.stride == 1 && inp == oup;

                    var layers = new List<nn.Module<Tensor, Tensor>>();
                    if (expand_ratio != 1) {
                        // pw
                        layers.Add(
                            Ops.Conv1dNormActivation(
                                inp,
                                hidden_dim,
                                kernel_size: 1,
                                norm_layer: norm_layer,
                                activation_layer: (inplace) => nn.ReLU6(inplace)));
                    }
                    layers.AddRange(new List<nn.Module<Tensor, Tensor>> {
                        // dw
                        Ops.Conv1dNormActivation(
                            hidden_dim,
                            hidden_dim,
                            stride: stride,
                            groups: hidden_dim,
                            norm_layer: norm_layer,
                            activation_layer: (inplace) => nn.ReLU6(inplace)),
                        // pw-linear
                        nn.Conv1d(hidden_dim, oup, 1, 1, 0, bias: false),
                        norm_layer(oup)
                    });
                    this.conv = nn.Sequential(layers);
                    this.out_channels = oup;
                    this._is_cn = stride > 1;
                    RegisterComponents();
                }

                public override Tensor forward(Tensor x)
                {
                    if (this.use_res_connect) {
                        return x + this.conv.call(x);
                    } else {
                        return this.conv.call(x);
                    }
                }
            }

            private readonly nn.Module<Tensor, Tensor> classifier;
            private readonly nn.Module<Tensor, Tensor> features;
            private readonly long last_channel;

            protected override void Dispose(bool disposing)
            {
                if (disposing) {
                    classifier.Dispose();
                    features.Dispose();
                }
                base.Dispose(disposing);
            }

            internal MobileNetV2(
                string name,
                long num_classes = 2,
                double width_mult = 1.0,
                long[][]? inverted_residual_setting = null,
                long round_nearest = 8,
                Func<long, long, long, long, Func<long, nn.Module<Tensor, Tensor>>, nn.Module<Tensor, Tensor>>? block = null,
                Func<long, nn.Module<Tensor, Tensor>>? norm_layer = null,
                double dropout = 0.2) : base(name)
            {
                if (block == null) {
                    block = (input_channel, output_channel, stride, t, norm_layer) => new InvertedResidual("InvertedResidual", input_channel, output_channel, stride, t, norm_layer);
                }

                if (norm_layer == null) {
                    norm_layer = (features) => nn.BatchNorm1d(features);
                }

                long input_channel = 32;
                long last_channel = 1280;
                
                

                if (inverted_residual_setting == null) {
                    inverted_residual_setting = new long[][] {
                        // t, c, n, s
                        new long[] { 6, 16, 1, 2 },
                        new long[] { 6, 24, 2, 2 },
                        new long[] { 6, 32, 3, 2 },
                        new long[] { 6, 64, 4, 2 }
                        // new long[] { 6, 96, 3, 1 },
                        // new long[] { 6, 160, 3, 2 },
                        // new long[] { 6, 320, 1, 1 }
                    };
                }

                // only check the first element, assuming user knows t,c,n,s are required
                if (inverted_residual_setting.Length == 0 || inverted_residual_setting[0].Length != 4) {
                    throw new ArgumentException($"inverted_residual_setting should be non-empty or a 4-element list, got {inverted_residual_setting}");
                }

                // building first layer
                input_channel = make_divisible.Make_divisible(input_channel * width_mult, round_nearest);
                this.last_channel = make_divisible.Make_divisible(last_channel * Math.Max(1.0, width_mult), round_nearest);
                var features = new List<nn.Module<Tensor, Tensor>> {
                    Ops.Conv1dNormActivation(1, input_channel, stride: 2, norm_layer: norm_layer, activation_layer: (inplace) => nn.ReLU6(inplace))
                };
                // building inverted residual blocks
                foreach (var x in inverted_residual_setting) {
                    long t = x[0];
                    long c = x[1];
                    long n = x[2];
                    long s = x[3];
                    var output_channel = make_divisible.Make_divisible(c * width_mult, round_nearest);
                    for (var i = 0; i < n; i++) {
                        var stride = i == 0 ? s : 1;
                        features.Add(block(input_channel, output_channel, stride, t, norm_layer));
                        input_channel = output_channel;
                    }
                }
                // building last several layers
                features.Add(
                    Ops.Conv1dNormActivation(
                        input_channel,
                        this.last_channel,
                        kernel_size: 1,
                        norm_layer: norm_layer,
                        activation_layer: (inplace) => nn.ReLU6(inplace)));
                // make it nn.Sequential
                this.features = nn.Sequential(features);

                // building classifier
                this.classifier = nn.Sequential(
                    nn.Dropout(p: dropout),
                    nn.Linear(this.last_channel, num_classes));

                RegisterComponents();

                // weight initialization
                foreach (var (_, m) in this.named_modules()) {
                    if (m is TorchSharp.Modules.Conv1d) {
                        var conv = (TorchSharp.Modules.Conv1d)m;
                        nn.init.kaiming_normal_(conv.weight, mode: nn.init.FanInOut.FanOut);
                        if (conv.bias is not null) {
                            nn.init.zeros_(conv.bias);
                        }
                    } else if (m is TorchSharp.Modules.BatchNorm1d) {
                        var norm = (TorchSharp.Modules.BatchNorm1d)m;
                        nn.init.ones_(norm.weight);
                        nn.init.zeros_(norm.bias);
                    } else if (m is TorchSharp.Modules.GroupNorm) {
                        var norm = (TorchSharp.Modules.GroupNorm)m;
                        nn.init.ones_(norm.weight);
                        nn.init.zeros_(norm.bias);
                    } else if (m is TorchSharp.Modules.Linear) {
                        var linear = (TorchSharp.Modules.Linear)m;
                        nn.init.normal_(linear.weight, 0, 0.01);
                        nn.init.zeros_(linear.bias);
                    }
                }
            }

            public override Tensor forward(Tensor x)
            {   
                x = TorchSharp.torch.reshape(x, -1, 1, 1000);
                x = this.features.call(x);
                // Cannot use "squeeze" as batch-size can be 1
                x = nn.functional.adaptive_avg_pool1d(x, 1);
                x = TorchSharp.torch.flatten(x, 1);
                x = this.classifier.call(x);
                return x;
            }
        }
//     }
// }