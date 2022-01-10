# Ray tracing

[Original repository](https://github.com/jevonlongdell/raytracing)

Implements [General Ray-Tacing Procedure](https://github.com/CSChisholm/raytracing/blob/master/Spencer-and-Murty_JOSA_52-6_1962.pdf)

## Installation

### Prerequisites
- numpy
- matplotlib
- scipy

### Download and install
```shell
git clone https://github.com/CSChisholm/raytracing.git
cd raytracing
python3 setup.py install
```

## Use
```python
from raytracing import raytracing as rt
```

See [Examples](https://github.com/CSChisholm/raytracing/tree/master/raytracing/Examples).

The example scripts can be imported and run using
```python
from raytracing import Examples as rtExamples
rtExamples.Imaging.main() 
```
