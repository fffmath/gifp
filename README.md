# GIFP

Code for the paper [“Generalized Implicit Factorization Problem"](https://eprint.iacr.org/2023/1562.pdf).

## Introduction

This is a Python implementation of GIFP based on [Joachim Vandersmissen's crypto-attacks](https://github.com/jvdsn/crypto-attacks) and [Mengce Zheng's Boneh_Durfee_Attack](https://github.com/MengceZheng/Boneh_Durfee_Attack).

## Requirements

- [SageMath](https://www.sagemath.org/) with Python 3.11.1. SageMath 9.8 is recommended.

  You can check your SageMath Python version using the following command:

  ```bash
  $ sage -python --version
  Python 3.11.1
  ```
Note: If your SageMath Python version is older than 3.11.1, some features in some scripts might not work.
## Usage
```bash
sage gifp.sage <modulus_bit_length> <alpha> <gamma> <beta1> <beta2> <m>
```
For example:
```bash
sage gifp.sage 200 0.1 0.7 0.1 0.15 4
```
If the roots are successfully found, it returns 1; otherwise, it returns 0. The corresponding time for computing LLL and computing Grobner basis can be found in the `gifp.log` file.

### Debug

You can enable debugging by setting `logging.basicConfig(filename='identifying_ideal_lattice.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')` in your code.

### Precautions
It is important to note that our paper introduces a new variable, 'w', to eliminate some 'z'. Introducing multiple variables may intuitively make it more challenging to satisfy the Grobner basis Heuristic. However, in practice, it is not necessary to satisfy the Grobner basis Heuristic to find the desired 'p' and 'q'. At the same time, if we abandon the introduction of 'w', the corresponding bound changes from $\gamma>4a(1-\sqrt{\alpha})$ to $\gamma>2a(2-\sqrt{\alpha})$. Even with only three variables in this case, we can still find 'p' and 'q' without satisfying the Grobner basis Heuristic.

And here's an example to illustrate the analysis above.

#### Example
Run the `example.sage`:
```bash
N1: 814072953237453703269626775081889594665601447233336309048661
N2: 1053039258722741964842067475826840037602528093495734686029849
The desired solution is $(x_0, y_0, z_0, w_0) = (-52202915, 633281, 809747, 1300454658952415958122805611909448306202465823887874467)$.
```
We can find the Groebner basis produced by the algorithm in the log:
```bash
DEBUG - Groebner basis x^2 + 104405830/633281*x*y + 2725144334497225/401044824961*y^2
DEBUG - Groebner basis x*y*z - 512797389907*x + 52202915/633281*y^2*z - 42271153812505*y
DEBUG - Groebner basis y^2*z^2 - 1025594779814*y*z + 262961163095431785468649
```
By calculating the third equation, we find that it is actually $(yz - 512797389907)^2$, thus $y_0z_0= 512797389907$.

Next, compute $z_0=q_2=\gcd(y_0z_0, N_2)=\gcd(512797389907, 1053039258722741964842067475826840037602528093495734686029849)= 809747$, hence $p_2=1300454658952415958122805611909448306202465823887874467$. $y_0=512797389907/809747=633281$.

Finally, we can compute from the first equation that $(x+52202915y/633281)^2=0$, thus $x_0=-52202915$.

### Author

You can find more information on [my personal website](https://www.fffmath.com/).

### License

This script is released under the MIT License. See the [LICENSE](LICENSE) file for details.