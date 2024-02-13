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

You can enable debugging by setting `logging.basicConfig(filename='gifp.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')` in your code.

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

Next, compute $z_0=q_2=\gcd(y_0z_0, N_2)=\gcd(512797389907, 1053039258722741964842067475826840037602528093495734686029849)= 809747$, hence $p_2=1300454658952415958122805611909448306202465823887874467$ and $y_0=512797389907/809747=633281$.

Finally, we can compute from the first equation that $(x+52202915y/633281)^2=0$, thus $x_0=-52202915$.

In fact, we can use this idea to solve the equations in GIFP with the addtional variable 'w'. We need to add $z*w-N_2$ to the polynomials set, which is used to compute grobner basis. The results are as follow:
```bash
DEBUG - Groebner basis x^3 + 949054662290877183643701285788447485524228569374337758441/1300454658952415958122805611909448306202465823887874467*x^2*w - 1420334021483155590500436908946912519349549371929066972100396907413475503310567648869804939900703343959343104/1300454658952415958122805611909448306202465823887874467*x*y*w + 899470549458874255507707186134809707166201975804710547998443302657499563269034401308056512214152546504854128863579/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*x*w^2 + 3351951982485649274893506249551461531869841455148098344430890360930441007518386744200468574541725856922507964546621512713438470702986642486608412251521024*y^3 - 4242542320356159613723265834226948229596340730191705146444381133524602217522250756879829992614500764308412951710072556854222740843066500782515326238488197595136/1300454658952415958122805611909448306202465823887874467*y^2*w + 1340593898074035779244258664229375579142604039668920409491792222011772083843752654270532390610853319333813450204404294352610256399531879122871090672283032242702254080/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*y*w^2 + 898238038084017563562338227646432387099377320303813307800488091356765552376408791836233455218847553585908501771721/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*w^3
DEBUG - Groebner basis x^2*y - 633281/1300454658952415958122805611909448306202465823887874467*x^2*w + 949054662290877183643701285788447485524228569374285555526/1300454658952415958122805611909448306202465823887874467*x*y*w - 601018285590228993735066793965393812080268992641916930889060806/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*x*w^2 - 2239744742177804210557442280568444278121645497234649534899989100963791871180160945380877493271607115776*y^3 + 2834829348730150494297645905695082363329911974473701304232603086128877736607962862085437499548906533697880064/1300454658952415958122805611909448306202465823887874467*y^2*w - 895773015334304179671600310669677746965728009302019031067149736696175450176300332988340217972126149493913907504695/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*y*w^2 - 600194732363352948804690825504678184747468788882499393858724579/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*w^3
DEBUG - Groebner basis x*y^2 - 1266562/1300454658952415958122805611909448306202465823887874467*x*y*w + 401044824961/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*x*w^2 + 1496577676626844588240573268701473812127674924007424*y^3 - 1894207960604897119413034154741166626129849741276803081821/1300454658952415958122805611909448306202465823887874467*y^2*w + 598547625909600858943938888583246930081868381363631260683838010/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*y*w^2 + 401044824961/1691182319991044503095767145359346461066521997444639731559007094215960447978904619750201603152003751150534089*w^3
DEBUG - Groebner basis z*w - 1053039258722741964842067475826840037602528093495734686029849
```
First, we compute $\gcd(f_i, f_j)$ as follow:
$$
\gcd(f_1,f_3)=x + 1496577676626844588240573268701473812127674924007424*y + w,\\
\gcd(f_1,f_2)=x + 1496577676626844588240573268701473812127674924007424*y + w,\\
\gcd(f_2,f_3)=1300454658952415958122805611909448306202465823887874467*x*y + 1946231412053562234409304066247558456395589794889005706640758220782555983621440415846687067372309088043008*y^2 - 633281*x*w - 946453752972972351727455674564628588911823637726457603677*y*w - 633281*w^2
$$
However, we find $\gcd(f_1,f_2)\vert \gcd(f_1,f_3)$, which yield 
$$
\gcd(f_2,f_3)=(1300454658952415958122805611909448306202465823887874467*y - 633281*w) * (x + 1496577676626844588240573268701473812127674924007424*y + w).
$$
Then we have $y_0z_0= 512797389907$ and the remaining steps are the same as above.
### Author

You can find more information on [my personal website](https://www.fffmath.com/).

### License

This script is released under the MIT License. See the [LICENSE](LICENSE) file for details.
