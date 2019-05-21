# System Idntification Library
系统辨识库
---

## Require
系统要求
---

+ Windows
+ Visual Studio 2017
+ [Eigen](https://github.com/eigenteam/eigen-git-mirror)
+ [polynomial](https://github.com/onathanventura/polynomial)

Lower version of VS may not open the .sln file but you can create your solution.

低版本VS可能无法打开工程，但是可以自己重新创建。

## M Sequence Generator
M序列生成器
---
M Sequence is a frequently used sequence in system identification to take the place of random.

M序列是常用于替代随机数发生器的用于系统辨识的生成序列。

This library contains two classes: a general sequence generator and a M sequence generator, where the later one is drived from the first.

本库中包含两个类：一个通用的序列生成器和一个M序列生成器，其中后者继承自前者。

For the general sequence generator, you need to specify the bit width, the coefficients of polynomial of the sequence and the initial state. Then you can get your sequence with get().

对于通用序列生成器，你需要指定位宽、序列多项式系数和初始状态。然后你就可以用get()获取序列了。

For M sequence generator, you only need to specify the bit width, then the coefficients and initial state will be auto settled.

对于M序列生成器，只需要指定位宽，随后系数和初始状态会自动设定。

## System Identification
系统辨识库
---
System identification library currently contains Model.h and LeastSquareEstimation.h.

系统辨识库目前包括Model.h和LeastSquareEstimation.h。

### Model.h
Here defines discrete domain model using ARX model in ArxModel class and continuous domain model in SModel class.

这里在ArxModel类中使用ARX模型定义了离散域模型，在SModel类中定义了连续域模型。

SModel class is simple. It only contains two Eigen::VectorXd representing the coefficients of numerator and denominator.

SModel类非常简单，它只包含了两个Eigen::VectorXd分别用来代表分子和分母的系数。

ArxModel contains a lot of parameters including order and delay time. Because this class can calculate the output with certain input sequence using the method update(double), I make its parameters protected to prevent from bing modified illegally form outside. It provides necessary interface to the parameters.

ArxModel包含了许多参数包括阶次、。由于这个类可以通过update(double)方法来计算某输入序列下的输出，我令它的参数为protected以避免被外部非法地修改。

Here also provides d2c() and c2d() to make a conversion between ArxModel and SModel. The calling form of these two functions is similar. Given a certain model, and you can also specify the sample time and convert method, the function will return the new model.

这里也提供了ArxModel与SModel之间的转换函数d2c() 和 c2d()。这两个函数的调用形式类似。给定一个模型，同时也可以额外指定采样时间和转换方法，函数会返回新的模型。

Convert methods are defined in enum DiscretizationMethod, which currently only support Zoh and Tustin. Zoh is impliemented with first order backward difference other then the true zero order holder. The polynomial library is used to simulate the fractional operation in the calculation process.

转换方法定义在枚举DiscretizationMethod中，目前仅支持Zoh和Tustin。Zoh使用一阶向后差分实现，而不是真正的零阶保持器。计算过程中使用了polynomial库来模拟分式运算。

### LeastSquareEstimation.h
Here provides least square algorithm and recursion least square algorithm class.

这里提供了最小二乘法和递推最小二乘法类。

leastSquare function estimate arx model with least square algorithm. Input data could be std::vector\<double> or double[]. Function will return an ArxModel object. In the calculation process we use Jacobi Svd to calculate the least square.

leastSquare使用最小二乘法估计arx模型。输入数据可以是std::vector\<double>或者double数组。计算过程中使用Jacobi SVD分解来计算最小二乘解。

RLS class provides Recursion Least Square algorithm with forgetting factor and damping factor. RLS has to record inputs and outputs, so make it class instead of function.

RLS类提供带有遗忘因子和阻尼因子的递推最小二乘法。RLS算法需要记录输入和输出，因此把它做成了类而不是函数。

Initialize the RLS object with the order of ArxModel and the delay time, then update RLS with system input and output value, and you can get estimated model with getModel() method. Besides, update() method will return true when the disturbance of coefficients is less than error threshold.

用Arx模型的阶数以及延迟环节参数初始化RLS对象，然后用系统的输入和输出更新模型，随后可以用getModel()方法获取拟合的模型。另外，getModel()方法会在误差小于误差限的时候返回true。