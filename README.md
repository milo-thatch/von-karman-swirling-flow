## von Kármán swirling flow

This Fortran program solves the second-order, nonlinear ODEs that govern the self-similar von Kármán swirling flow, an interesting textbook exercise.

$$
2F + H^\prime = 0,
$$

$$
F^2 - G^2 + F^\prime H = F^{\prime\prime},
$$

$$
2FG + G^\prime H = G^{\prime\prime},
$$

To do so, the code splits the governing equations into a system of first-order ODEs with the auxiliary variables $A=F^\prime$ and $B=G^\prime$. A block-elimination (Thomas) algorithm is then used to solve the system iteratively, the nonlinearity being treated with the Newton-Raphson method.

The main assumptions:
<ul>
  <li>steady state, incompressible flow</li>
  <li>infinite plane disk</li>
  <li>azimuthal symmetry</li>
  <li>no rotation at infinity</li>
</ul>

To run, first install gfortran and LAPACK by running the commands
```console
sudo apt-get install gfortran
sudo apt-get install liblapack-dev
```
then navigate to the code's directory and allow your machine to run the bash script provided
```
chmod +x deploy.sh
```
and run by typing
```
./deploy.sh
```
The input.in file is used to specify the number of points and domain width.

## References
Aref, H., & Balachandar, S. (2018). *A first course in computational fluid dynamics*. Cambridge University Press. ISBN: 9781107178519.

Cebeci, T. (2002). *Convective heat transfer*. Horizons Pub. LCCN: 2002068512.



