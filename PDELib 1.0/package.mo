package PDE

  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>
This package contains necessary blocks for solving partial differential equations with two numerical methods: Method of Lines and Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));

  package UsersGuide "User큦 Guide"
    annotation (DocumentationClass=true, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Users Guide of the PDE Library</font></h3>
<p>
PDE Package provides necessary blocks that can be used together with standard Modelica block to implement partial differential equations. <br>
The numerical methods for solving PDEs that are implemented in Package are the Method of Lines (MOL) and the Finite Volume Method (FVM). <br>
In the following a short introduction to both methods is done. For detailled information see the documentation of the single blocks <br>
and Examples of the corresponding methods.
</p>

<ul>
<li> <a href=\"Modelica://PDE.UsersGuide.MethodOfLines\">Method of Lines</a> </li>
<li> <a href=\"Modelica://PDE.UsersGuide.FiniteVolumeMethod\">Finite Volume Method</a> </li>
<li> <a href=\"Modelica://PDE.UsersGuide.FluxLimiter\">Flux Limiter</a> </li>
</ul>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));

    package MethodOfLines
      annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Method of Lines</font></h3>
<p>
The first step in constructing PDE is to drag the <b>WorldModel</b>
block in the diagram. The WorldModel block contains general <br>
information about the problem, such as the number of grid points the
user wish to use to solve the problem. These information are then <br>
\"propagated\" to all the blocks in the diagram that need them.
The heart of the Method of Lines Package are the integrator and <br>
derivative blocks. Once the user know the problem it is easy to
implement it in Dymola. Typically the problem consists of a <br>
PDE and initial and boundary conditions. The initial condition can
be passed to the <b>IC</b> input of the integrator block. For the <br>
boundary conditions we must additionally specify some parameters in
the <b>WorldModel</b> block: <br><br>

<ul>
<li> <b>bcl</b> and <b>bcr</b> specify whether there is a boundary condition at
     the left and at the right respectively (0: no, 1: yes). <br><br> </li>
<li> <b>vb</b> and <b>ve</b> specify the first and last unknown variable
     respectively. The same applies to <b>icb</b> and <b>ice</b>, which specify the <br>
     first and last variable for which an initial condition must be
     specified.</li>
</ul>

<br>
The information in the second point is redundant and could be
deduced from the first point. However, for better understanding and <br>
use this choice was taken. In the WorldModel block the number
of grid points can be specified in the <b>n</b> parameter. Per <br>
default n = 10. Let now turn our attention to the integrator
block. This block implements PDEs of the form

<p>
<img align=middle src=\"..\\Images\\ref1.png\">
</p>

where <b>u</b> is the unknown function and <b>R</b> is the right
part of the equation that can contain space derivatives, constant <br>
values and so on (See Examples). This means that if we do not have
the PDE in this form we must first convert it to this form before <br>
passing it to the integrator. Once this step is achieved we can
build the right part of the equation and at the end pass it to the <br>
<b>R</b> input of the integrator block.
Now say that one of the component of the right part of the equation <br>
is a derivative of unknown function with respect to space. In this
case we can use the derivative blocks that the Package provide. For <br>
example for the first-order space derivative use the <b>u_x</b>
block. It is important to note here that in the case of the boundary <br>
condition of type

<p>
<img align=middle src=\"..\\Images\\ref2.png\">
</p>

say at the left part of the domain, we must set the bcl
parameter in the corresponding block to -1 which tells that we are <br>
in front of the symmetry boundary condition at the left part
of the domain. The same applies for the right part of the domain, in <br>
which case we must set the parameter bcr to -1. Another
important point is the accuracy of the derivative computation. For <br>
the computation of the first order space derivative for example, the
second, fourth and sixth order central difference approximations <br>
were provided. The user can choose which one to use by specifying
the value of the <b>u_x</b> parameter in the <b>WorldModel</b> block. <br>

</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end MethodOfLines;

    package FiniteVolumeMethod
      annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Finite Volume Method</font></h3>
<p>
Before going into details of FVM Package we need to explain some
theory. In the following a short introduction to the finite volume <br>
method is made. This short introduction is based on book of Randall
Leveque, \"Finite Volume Methods for Hyperbolic Problems\".
<ul>
<li> <a href=\"Modelica://PDE.UsersGuide.FiniteVolumeMethod.Introduction\">Introduction</a> </li>
<li> <a href=\"Modelica://PDE.UsersGuide.FiniteVolumeMethod.UnstableMethod\">Unstable Method</a> </li>
<li> <a href=\"Modelica://PDE.UsersGuide.FiniteVolumeMethod.LaxFriedrichMethod\">Lax-Friedrichs Method</a> </li>
<li> <a href=\"Modelica://PDE.UsersGuide.FiniteVolumeMethod.Implementation\">Implementation</a> </li>
</ul>


</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      package Introduction
        annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Introduction</font></h3>
<p>

  In one dimension, the finite volume method consists in subdividing the spatial domain
  into intervals, \"finite volumes\" (or cells), and approximate the integral of the function <b>q</b> over each <br>
  of these volumes at each time step.  <br>
  Denote the i-th finite volume by

<p>
<img align=middle src=\"..\\Images\\ref3.png\">
</p>

  Then the approximation to the average of <b>q</b> in the cell C<sub>i</sub> at time
  t, which we denote with Q<sub>i</sub><sup>t</sup>, is

<p>
<img align=middle src=\"..\\Images\\ref4.png\">
</p>

  Remains the question of how to find this approximation. If we
  think about conservation law, we note that the average within the
  cell can only changes due to the fluxes at the boundaries <br>
  (if we assume that no source or sink is present in the cell).
  The integral form of conservation law is

<p>
<img align=middle src=\"..\\Images\\ref5.png\">
</p>

  If we integrate this expression in time from t to t+Delta_t,
  we obtain

<p>
<img align=middle src=\"..\\Images\\ref6.png\">
</p>

  and dividing by Delta_x we reach the form

<p>
<img align=middle src=\"..\\Images\\ref7.png\">
</p>

  which gives us an explicit time marching algorithm. This is more
  clearly seen if we rewrite the expression using the notation we
  introduced above:

<p>
<img align=middle src=\"..\\Images\\ref8.png\">
</p>

  where F<sub>i-1/2</sub><sup>t</sup> approximates the average flux along the
  interface x<sub>i-1/2</sub>:

<p>
<img align=middle src=\"..\\Images\\ref9.png\">
</p>

  As can be seen from the equation (??), in order to find the average
  at the next time step we need to find the fluxes at the
  interfaces. The flux at the interface x<sub>i-1/2</sub> for example, <br>
  depends on q(x<sub>i-1/2</sub>, t), which changes with time along the interface and for
  which we do not know the analytical solution. For this reason we need to find
  some approximation to <br>
  these fluxes in order to calculate the averages at the next
  time step. Let us now see some simple flux approximations.
  Examples of fluxes: <br><br>

  <b>Advection equation</b>
  <br> <br>

  Consider the advection equation q<sub>t</sub> + uq<sub>x</sub> = 0,
  where <b>u</b> is the fluid velocity. We have seen in the
  previous chapters, that the flux of <br>
  the contaminant at some point x, at some time t, could be written as uq(x, t).
  Consider now the flux through the interface x<sub>i-1/2</sub>. <br>

  Inserting it into the average update rule, we obtain the finite volume method for
  the advection equation:

<p>
<img align=middle src=\"..\\Images\\ref10.png\">
</p>

<br><br>

  <b>Diffusion equation</b>
  <br><br>

  In the advection equation, the flux depends on <b>q: f(q) =
  uq</b>. The flux in the diffusion equation depends on
  the derivative of <b>q</b>:

<p>
<img align=middle src=\"..\\Images\\ref11.png\">
</p>

  where <b>beta</b> is the conductivity. If <b>beta</b> is space dependent
  then the flux will depend on space too (<b>f(x, q<sub>x</sub>) = -beta(x)
  q<sub>x</sub></b>). In the following we will assume for simplicity <br>
  that <b>beta</b> is constant.
  Now remains the question of how to approximate numerically the
  diffusion flux. One possibility were:

<p>
<img align=middle src=\"..\\Images\\ref12.png\">
</p>

  By inserting this flux approximation into the average update rule
  (??), we obtain:

<p>
<img align=middle src=\"..\\Images\\ref13.png\">
</p>

  It is interesting to note, that after some algebraic
  manipulations, we can write the average update rule
  in the form

<p>
<img align=middle src=\"..\\Images\\ref14.png\">
</p>

  which is equivalent to the finite difference discretization
  of the conservation law equation <b>q<sub>t</sub> + f(q)<sub>x</sub> = 0</b>.
  As said in (LeVeque): Many methods can be equally well viewed as <br>
  finite difference approximations to this equation or as
  finite volume methods.
  Another form of the average update rule is

<p>
<img align=middle src=\"..\\Images\\ref15.png\">
</p>

  which give us an ODE for each average cell. This form is more
  suitable for the implementation in Dymola and all Finite Volume
  Method blocks are based on this form of update rule.

</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Introduction;

      package UnstableMethod
        annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Unstable Method</font></h3>
<p>

The unstable flux just take the arithmetic average of the
fluxes based on either side of the interface x<sub>i-1/2</sub>:

<p>
<img align=middle src=\"..\\Images\\ref16.png\">
</p>

By using this flux, the average update rule becomes:

<p>
<img align=middle src=\"..\\Images\\ref17.png\">
</p>

This method is generally unstable.

</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end UnstableMethod;

      package LaxFriedrichMethod
        annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Lax-Friedrichs Method</font></h3>
<p>

The Lax-Friedrichs flux is defined as:

<p>
<img align=middle src=\"..\\Images\\ref18.png\">
</p>

inserting it into the average update rule, we obtain the
Lax-Friedrichs method:

<p>
<img align=middle src=\"..\\Images\\ref19.png\">
</p>

If we take a closer look at the Lax-Friedrichs flux, we notice
that it is similar to the unstable flux, but with the addition
of some correction term. This correction term looks like
the diffusion flux (??) with

<p>
<img align=middle src=\"..\\Images\\ref20.png\">
</p>

The Lax-Friedrichs flux can thus be interpreted as the unstable
flux plus a numerical diffusion. This numerical diffusion damps
the instabilities that arise in the unstable method, however,
it damps it too much. <br>
Later we will see another method, the
Lax-Wendroff method, that add less diffusion.

</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end LaxFriedrichMethod;

      package Implementation
        annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Implementation</font></h3>
<p>

  As already seen, we obtain by starting from the integral
  conservation law the cell average update rule

<p>
<img align=middle src=\"..\\Images\\ref21.png\">
</p>

  This ODE is implemented in <b>FVMIntegrator</b> block.
  The use of this block is similar to the use of the integrator
  block in the MOL package. <br>
  We must specify the initial and boundary conditions
  and give them to the inputs <b>IC</b> and <b>gcl</b>, <b>gcr</b>
  respectively. The only difference here is that <br>
  instead of writing special formulas for the boundary conditions like in the MOL
  package, we extend the domain with additional cells, called <b>ghost
  cells</b>, <br>
  and fill them with values. This way, the first and the last
  cell of the domain for example, can use the same formula as all
  other cells in the domain. <br>
  Finally, we must pass the \\emph{flux}
  vector <b>F</b> to the <b>F</b> input of the integrator, which can
  now use the cell average update rule with the information <br>
  just provided and compute the averages for the next time step.
  See Examples implemented in the package for better understanding.

</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Implementation;
    end FiniteVolumeMethod;

    package FluxLimiter
      annotation (Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Flux Limiter</font></h3>
<p>

Given a hyperbolic system of <b>m</b> equations we want to solve it
with the flux limiter method. We can specify the parameter <b>m</b>
as well as the number of cell averages <br>
<b>n</b> that we wish to use in the <b>WorldModel</b> block.
The first block we need is the integrator block that will
give the <b>mxn+gcl+gcr</b> average matrix <b>Q</b> as output <br>
with which we can start the construction. Having matrix
<b>Q</b> we can compute the jumps

<p>
<img align=middle src=\"..\\Images\\ref22.png\">
</p>

at each interface i = 1,...,n+gcl+gcr-1. The <b>deltaQ</b> block
achieves this task. The next step is to solve the Riemann problem

<p>
<img align=middle src=\"..\\Images\\ref23.png\">
</p>

This can be accomplished by passing <b>DeltaQ</b> and <b>mxm</b>
eigenvalue matrix <b>R</b> to the <b>Riemann</b> block which will
give us the <br>
<b>mxn+gcl+gcr-1</b> matrix <b>alpha</b> as output. It is
important to note that the eigenvalue matrix <b>R</b> as well as the
eigenvalues <b>lambda<sub>i</sub></b> must be provided <br>
by the user. In the case of a constant coefficient linear hyperbolic system these do not
change with time. The next step is to use the <b>alpha</b> matrix just
obtained together <br>
with the <b>i-th</b> eigenvalue <b>lambda<sub>i</sub></b> to
calculate <b>theta<sub>i</sub></b>. The <b>theta<sub>i</sub></b> matrix has the <b>theta<sub>i</sub></b>
values in the <b>i-th</b> row and zeros elsewhere. <br>
Once <b>theta<sub>i</sub></b> matrix is computed we can use for instance <b>Beam-Warming</b> block to
compute <b>phi(theta)</b>, which is just <b>theta</b> in the Beam-Warming
method. <br>
Beam-Warming together with many other methods is implemented
in <b>FluxSolver</b> block. The user can choose which method to use
in the <b>WorldModel</b> block <br>
by giving the corresponding value to
the <b>fls</b> variable. Here is a list of methods with their
corresponding values:

<br> <br>
<pre>
<ul>
<li> <b><font color=\"#33CCFF\" size=3>Upwind method:</font></b>          fls = 1 </li>
<li> <b><font color=\"#33CCFF\" size=3>Lax-Wendroff method:</font></b>    fls = 2 </li>
<li> <b><font color=\"#33CCFF\" size=3>Beam-Warming method:</font></b>    fls = 3 </li>
<li> <b><font color=\"#33CCFF\" size=3>Fromm method:</font></b>           fls = 4 </li>
<li> <b><font color=\"#33CCFF\" size=3>van Leer method:</font></b>        fls = 5 </li>
</ul>
</pre>

The next step is to pass the <b>alpha</b> and <b>R</b> matrix to the
<b>pWave</b> block, that will calculate the <b>p-th</b> Wave matrix

<p>
<img align=middle src=\"..\\Images\\ref24.png\">
</p>

and give the <b>mxn+1</b> matrix as output. Each column <b>p</b> in
this matrix contains the wave <b>W<sub>i-1/2</sub><sup>p</sup></b>. By using the same
block but giving this time <br>
<b>widetilde{alpha}</b> as input instead of
<b>\\alpha</b> we can compute the limited wave <b>widetilde{W}</b>

<p>
<img align=middle src=\"..\\Images\\ref25.png\">
</p>

The limited <b>widetilde{alpha}</b> can be computed by using the block
<b>LimitedAlpha</b> which need <b>phi(theta)</b> and <b>alpha</b> matrices
as input. With the eigenvalues <br>
<b>lambda<sub>p</sub></b> and the waves matrices
<b>W</b>, <b>widetilde{W}</b> we can compute the fluxes and fluctuations.
This is done by blocks <b>FluxLimited</b> and <b>Fluctuation</b>. <br>
The <b>FluxLimited</b> block computes

<p>
<img align=middle src=\"..\\Images\\ref26.png\">
</p>

If we have m equations in the system then we must use m
<b>FluxLimited</b> block to compute each term of the sum

<p>
<img align=middle src=\"..\\Images\\ref27.png\">
</p>

and then sum all the outputs of the blocks and pass the result to the
<b>Flux</b> input of the integrator block. The same applies to the
<b>Fluctuation</b> block, only that here we must sum the <b>+</b> outputs <br>
of the blocks and <b>-</b> outputs of the blocks separately
and at the end pass them to the <b>+</b> and <b>-</b> inputs of the
integrator, respectively. Finally we specify the initial <br>
and boundary conditions and connect them to the <b>IC</b> and <b>gcl</b>,
<b>gcr</b> inputs, respectively.


</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end FluxLimiter;
  end UsersGuide;

  package Icons
    partial block BlockIcon
      annotation (Icon(Rectangle(extent=[-100,100; 100,-100], style(
              color=3,
              rgbcolor={0,0,255},
              thickness=2,
              fillColor=6,
              rgbfillColor={255,255,0}))));
    end BlockIcon;

    partial block BlockIcon1
      annotation (Icon(Bitmap(
            extent=[-100,133; 100,-133],
            style(thickness=2),
            name="Images/int1.png")));
    end BlockIcon1;

    partial block BlockIcon2
      annotation (Diagram, Icon(Bitmap(extent=[-100,100; 100,-100], name=
                "Images/d_dx1.png")));
    end BlockIcon2;

    partial block BlockIcon3
      annotation (Icon(Bitmap(extent=[-100,102; 100,-102], name=
                "Images/d2_dx21.png")));
    end BlockIcon3;

    partial block BlockIcon4
      annotation (Icon(Rectangle(extent=[-100,100; 100,-100], style(
              color=3,
              rgbcolor={0,0,255},
              thickness=4,
              fillColor=1,
              rgbfillColor={255,0,0}))));
    end BlockIcon4;

    partial block BlockIcon5
      annotation (Icon(Rectangle(extent=[-100,100; 100,-100], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=68,
              rgbfillColor={170,213,255}))));
    end BlockIcon5;
    annotation (Documentation(info="<html>
<p>
This package contains definitions for the graphical layout of
components which may be used in different libraries.
The icons can be utilized by inheriting them in the desired class
using \"extends\".
</p>
</html>"));
  end Icons;

  package World
    model worldModel
      extends PDE.Icons.BlockIcon;

    parameter Integer n = 10 "Number of grid points or cells";
    parameter Integer qss=1
        "First-order space derivative (qss 1: second-order, 2: fourth-order, 3: sixth-order)";
    parameter Integer u_xx=1
        "Second-order space derivative(u_xx Method 1: second-order, 2: fourth-order)";

    parameter Integer gcl = 2
        "|Finite Volume Method| Number of left ghost cells";
    parameter Integer gcr = 2
        "|Finite Volume Method| Number of right ghost cells";
    parameter Real deltax = 1/n "|Finite Volume Method| Space interval";
    parameter Real deltat = 0 "|Finite Volume Method| Time interval";

    parameter Integer m = 0
        "|Finite Volume Method| Flux Limiter: Number of unknowns in the system of equations";

    parameter Integer fls=1
        "|Finite Volume Method| Flux Limiter Method 1: Upwind, 2: Lax Wendroff, 3: Beam Warming, 4: Fromm, 5: van Leer";

      annotation (defaultComponentName="worldModel1", defaultComponentPrefixes="inner", Diagram,
            Icon(Text(
            extent=[-64,52; 70,-40],
            style(color=3, rgbcolor={0,0,255}),
            string="%name")),
        Documentation(info="<html>
<p>
The <b>worldModel</b> block provides information that are propagated to other blocks. If for example the user <br>
change the number of grid points, <i>n</i>, this change will propagate to all blocks in the model that need this
information. <br>
Other information that are propagated are listed in the table below. It is for example also possible to <br>
to chose the order of the derivation.
</p>
<p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end worldModel;
    annotation (Documentation(info="<html>
<p>
This package contains the worldModel block that provides global informations for other blocks.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
  end World;

  package MOL
    package SpaceDerivative

      annotation (Documentation(info="<html>
<p>
The SpaceDerivative package contains blocks for the computation of the derivatives with respect to space. <br>
The derivative computation is based on Newton-Gregory backward polynomials
</p>
<img align=middle src=\"..\\Images\\f1.png\">
<p>
where x<sub>n</sub> is the last grid point, s = (x-x<sub>n</sub>)/h and h = x<sub>i+1</sub> - x<sub>i</sub> (for i = 0,..., n-1). <br>
By using this polynomial we can compute the first and second space derivatives, that are implemented in the <br>
respective blocks. To compute for example the first derivative we must do
</p>
<img align=middle src=\"..\\Images\\f2.png\">
<p>
where h = x<sub>i+1</sub> - x<sub>i</sub> (for i = 0,..., n-1). We assume here that the grid points are spaced equally.
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      package Derivatives
        block u_x
          extends PDE.Icons.BlockIcon2;

          outer PDE.World.worldModel worldModel1;
          parameter Integer method = worldModel1.qss;
          inner parameter Integer n = worldModel1.n;

          inner parameter Integer bcl = 0
            "|Boundary Conditions| Type of the boundary condition at the left (-1:symmtery; 0: none)";
          inner parameter Integer bcr = 0
            "|Boundary Conditions| Type of the boundary condition at the right (-1:symmtery; 0: none)";

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 122,10]);
          annotation (Icon,           Diagram,
            Documentation(info="<html>
<p>
The <b>u_x</b> block computes the first-order space derivative. By using the Newton-Gregory backward polynomial we <br>
obtain
</p>

<img align=middle src=\"..\\Images\\f2.png\">

<p>
If we wish a second-order central difference approximation, we need to fit the polynomial through the three points <br>
x<sub>i-1</sub>, x<sub>i</sub>, x<sub>i+1</sub>. This means to write the polynomial for example around the point x<sub>i+1</sub>,
drop the higher-oder terms and set s = -1 to obtain
</p>

<img align=middle src=\"..\\Images\\f9.png\">

<p>
where h = x<sub>i+1</sub> - x<sub>i</sub> (for i = 0,..., n-1). We assume here that the grid points are spaced equally. <br>
For the boundary point x<sub>1</sub> we use biased formula and obtain
</p>

<img align=middle src=\"..\\Images\\f10.png\">

<p>
By using the same idea we obtain a biased formula for the boundary point x<sub>n</sub>
</p>

<img align=middle src=\"..\\Images\\f11.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          PDE.MOL.SpaceDerivative.SDInterfaces.u_xCD4B4 der5_1 if
                                                  method == 1
                                     annotation (extent=[-20,-10; 0,10]);
          annotation (extent=[-20,-10; 0,10]);
          PDE.MOL.SpaceDerivative.SDInterfaces.u_xCD2B2 der1 if
                                                method == 2
                                  annotation (extent=[-20,60; 0,80]);
          PDE.MOL.SpaceDerivative.SDInterfaces.u_xCD6B6 cD6B6_1 if
                                                   method == 3
            annotation (extent=[-20,-80; 0,-60]);
        equation
          connect(u, der5_1.u) annotation (points=[-120,0; -22,0],
              style(color=74, rgbcolor={0,0,127}));
          connect(der5_1.y, y) annotation (points=[1,0; 111,0],               style(
                color=74, rgbcolor={0,0,127}));
          connect(u, der1.u)
            annotation (points=[-120,0; -72,0; -72,70; -22,70],
                                                style(color=74, rgbcolor={0,0,127}));
          connect(u, cD6B6_1.u) annotation (points=[-120,0; -72,0; -72,-70; -22,
                -70], style(color=74, rgbcolor={0,0,127}));
          connect(cD6B6_1.y, y) annotation (points=[1,-70; 52,-70; 52,0; 111,0],
              style(color=74, rgbcolor={0,0,127}));
          connect(der1.y, y) annotation (points=[1,70; 52,70; 52,0; 111,0],
              style(color=74, rgbcolor={0,0,127}));
        end u_x;

        block u_xx
          extends Icons.BlockIcon3;
        outer PDE.World.worldModel worldModel1;
        inner parameter Integer n = worldModel1.n;
        parameter Integer u_xx = worldModel1.u_xx;

          inner parameter Integer bcl = 0
            "|Boundary Conditions| Type of the boundary condition at the left (-1:symmtery; 0: none)";
          inner parameter Integer bcr = 0
            "|Boundary Conditions| Type of the boundary condition at the right (-1:symmtery; 0: none)";

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon,
            Documentation(info="<html>
<p>
The <b>u_xx</b> block computes the second-order space derivative. By using the Newton-Gregory backward polynomial we <br>
obtain
</p>

<img align=middle src=\"..\\Images\\f3.png\">

<p>
where h = x<sub>i+1</sub> - x<sub>i</sub> (for i = 0,..., n-1). We assume here that the grid points are spaced equally. <br>
If we wish a second-order central difference approximation, we need to fit the polynomial through the three points <br>
x<sub>i-1</sub>, x<sub>i</sub>, x<sub>i+1</sub>. This means to write the polynomial for example around the point x<sub>i+1</sub> and drop the higher-oder terms to obtain
</p>

<img align=middle src=\"..\\Images\\f4.png\">

<p>
and finally, to evaluate the second-order space derivative around the point x<sub>i</sub> we need to set s = -1 to get
</p>

<img align=middle src=\"..\\Images\\f5.png\">

<p>
The second-order central difference scheme is implemented in <b>u_xxCD2B2</b> block. <br>
By following the same approach we can compute the fourth-order central difference scheme. This time we need more
terms in the polynomial
</p>

<img align=middle src=\"..\\Images\\f6.png\">

<p>
and so we obtain
</p>

<img align=middle src=\"..\\Images\\f7.png\">

<p>
for the boundary points we use a biased formula and we obtain
</p>

<img align=middle src=\"..\\Images\\f8.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          PDE.MOL.SpaceDerivative.SDInterfaces.u_xxCD2B2 u_xxCD2B2_1 if
                                   u_xx == 1 annotation (extent=[-20,16; 12,40]);
          PDE.MOL.SpaceDerivative.SDInterfaces.u_xxCD4B4 derivatorSecond if
                                             u_xx == 2
            annotation (extent=[-20,-40; 12,-16]);
        equation
          connect(u, u_xxCD2B2_1.u) annotation (points=[-120,0; -72,0; -72,28;
                -23.2,28], style(color=74, rgbcolor={0,0,127}));
          connect(u, derivatorSecond.u) annotation (points=[-120,0; -72,0; -72,
                -28; -23.2,-28], style(color=74, rgbcolor={0,0,127}));
          connect(u_xxCD2B2_1.y, y) annotation (points=[13.6,28; 60,28; 60,0;
                110,0],
                    style(color=74, rgbcolor={0,0,127}));
          connect(derivatorSecond.y, y) annotation (points=[13.6,-28; 60,-28;
                60,0; 110,0],
                           style(color=74, rgbcolor={0,0,127}));
        end u_xx;
        annotation (Documentation(info="<html>
<p>
The Derivatives package contains blocks for the computation of the derivatives with respect to space. <br>
The derivative computation is based on Newton-Gregory backward polynomials
</p>
<img align=middle src=\"..\\Images\\f1.png\">
<p>
where x<sub>n</sub> is the last grid point, s = (x-x<sub>n</sub>)/h and h = x<sub>i+1</sub> - x<sub>i</sub> (for i = 0,..., n-1). <br>
By using this polynomial we can compute the first and second space derivatives, that are implemented in the <br>
respective blocks. To compute for example the first derivative we must do
</p>
<img align=middle src=\"..\\Images\\f2.png\">

<p>
where h = x<sub>i+1</sub> - x<sub>i</sub> (for i = 0,..., n-1). We assume here that the grid points are spaced equally.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Derivatives;

      package SDInterfaces
        block u_xCD2B2
          extends PDE.Icons.BlockIcon;

        outer Integer n;
        outer Integer bcl;
        outer Integer bcr;
        parameter Real deltax = 1/(n-1);
        // parameter Integer bcl = worldModel1.bcl;
        // parameter Integer bcr = worldModel1.bcr;

        equation
          if bcl == -1 then
            y[1] = 0;
          else
            y[1] = (-u[3] + 4*u[2] - 3*u[1])/(2*deltax);
          end if;

          for i in 2:n-1 loop
            y[i] = (u[i+1] - u[i-1])/(2*deltax);
          end for;

          if bcr == -1 then
            y[n] = 0;
          else
            y[n] = (3*u[n] - 4*u[n-1] + u[n-2])/(2*deltax);
          end if;

        public
          Modelica.Blocks.Interfaces.RealInput u[n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-70,36; 74,-34],
                style(color=3, rgbcolor={0,0,255}),
                string="CD2B2")),
            Documentation(info="<html>
<p>
Computes the second-order central difference approximation to the first-order space derivative with second-order biased approximation for boundary points. <br>
By using the Newton-Gregory backward polynomial we obtain<br>
</p>

<img align=middle src=\"..\\Images\\f2.png\">

<p>
For second-order central difference approximation, we need to fit the polynomial through the three points <br>
x<sub>i-1</sub>, x<sub>i</sub>, x<sub>i+1</sub>. This means to write the polynomial for example around the point x<sub>i+1</sub>,
drop the higher-oder terms and set s = -1 to obtain
</p>

<img align=middle src=\"..\\Images\\f9.png\">


<p>
For the boundary point x<sub>1</sub> we use biased formula and obtain
</p>

<img align=middle src=\"..\\Images\\f10.png\">

<p>
By using the same idea we obtain a biased formula for the boundary point x<sub>n</sub>
</p>

<img align=middle src=\"..\\Images\\f11.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end u_xCD2B2;

        block u_xCD4B4
          extends PDE.Icons.BlockIcon;

        outer Integer n;
        outer Integer bcl;
        outer Integer bcr;
        parameter Real deltax = 1/(n-1);
        // parameter Integer bcl = worldModel1.bcl;
        // parameter Integer bcr = worldModel1.bcr;

          Modelica.Blocks.Interfaces.RealInput u[n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[n]
            annotation (extent=[100,-10; 120,10]);
        equation
          if bcl == -1 then
            y[1] = 0;
            y[2] = (-u[4] + 8*u[3] + u[2] - 8*u[1])/(12*deltax);
          else
            y[1] = (-3*u[5] + 16*u[4] - 36*u[3] + 48*u[2] - 25*u[1])/(12*deltax);
            y[2] = (u[5] - 6*u[4] + 18*u[3] - 10*u[2] - 3*u[1])/(12*deltax);
          end if;

          for i in 3:n-2 loop
            y[i] = (-u[i+2] + 8*u[i+1] - 8*u[i-1] + u[i-2])/(12*deltax);
          end for;

          if bcr == -1 then
            y[n-1] = (8*u[n] - u[n-1] - 8*u[n-2] + u[n-3])/(12*deltax);
            y[n] = 0;
          else
            y[n-1] = (3*u[n] + 10*u[n-1] - 18*u[n-2] + 6*u[n-3] - u[n-4])/(12*deltax);
            y[n] = (25*u[n] - 48*u[n-1] + 36*u[n-2] - 16*u[n-3] + 3*u[n-4])/(12*deltax);
          end if;

          annotation (Icon(Text(
                extent=[-74,64; 74,-66],
                style(color=3, rgbcolor={0,0,255}),
                string="CD4B4")),Diagram,
            Documentation(info="<html>
<p>
Computes the fourth-order central difference approximation to the first-order space derivative with fourth-order biased approximation for boundary points.
</p>

<img align=middle src=\"..\\Images\\f12.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end u_xCD4B4;

        block u_xCD6B6
          extends Icons.BlockIcon;

        outer Integer n;
        outer Integer bcl;
        outer Integer bcr;
        parameter Real deltax = 1/(n-1);

        equation
          if bcl == -1 then
            y[1] = 0;
            y[2] = (u[5] - 9*u[4] + 44*u[3] + 9*u[2] - 45*u[1])/(60*deltax);
            y[3] = (u[6] - 9*u[5] + 45*u[4] - 46*u[2] + 9*u[1])/(60*deltax);
          else
            y[1] = (265*u[7] - 478*u[6] + 50*u[5] + 400*u[4] - 450*u[3] + 360*u[2] - 147*u[1])/(60*deltax);
            y[2] = (-3463*u[7] + 13845*u[6] - 20740*u[5] + 13760*u[4] - 3315*u[3] - 77*u[2] - 10*u[1])/(60*deltax);
            y[3] = (-u[7] + 8*u[6] - 30*u[5] + 80*u[4] - 35*u[3] - 24*u[2] + 2*u[1])/(60*deltax);
          end if;

          for i in 4:n-3 loop
            y[i] = (u[i+3] - 9*u[i+2] + 45*u[i+1] - 45*u[i-1] + 9*u[i-2] - u[i-3])/(60*deltax);
          end for;

          if bcr == -1 then
            y[n-2] = (-9*u[n] + 46*u[n-1] - 45*u[n-3] + 9*u[n-4] - u[n-5])/(60*deltax);
            y[n-1] = (45*u[n] - 9*u[n-1] - 44*u[n-2] + 9*u[n-3] - u[n-4])/(60*deltax);
            y[n] = 0;
          else
            y[n-2] = (-2*u[n] + 24*u[n-1] + 35*u[n-2] - 80*u[n-3] + 30*u[n-4] - 8*u[n-5] + u[n-6])/(60*deltax);
            y[n-1] = (10*u[n] + 77*u[n-1] - 150*u[n-2] + 100*u[n-3] - 50*u[n-4] + 15*u[n-5] - 2*u[n-6])/(60*deltax);
            y[n] = (147*u[n] - 360*u[n-1] + 450*u[n-2] - 400*u[n-3] + 225*u[n-4] - 72*u[n-5] + 10*u[n-6])/(60*deltax);
          end if;

        public
          Modelica.Blocks.Interfaces.RealInput u[n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Text(
                extent=[-76,54; 74,-54],
                style(color=3, rgbcolor={0,0,255}),
                string="CD6B6")), Documentation(info="<html>
<p>
Computes the sixth-order central difference approximation to the first-order space derivative with sixth-order biased approximation for boundary points.
</p>

<img align=middle src=\"..\\Images\\f13.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end u_xCD6B6;

        block u_xxCD2B2
          extends Icons.BlockIcon;

        outer Integer n;
        outer Integer bcl;
        outer Integer bcr;
        parameter Real deltax = 1/(n-1);
          Modelica.Blocks.Interfaces.RealInput u[n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[n]
            annotation (extent=[100,-10; 120,10]);
        equation
          //bcl = -1 means: symmetry condition. If we do not have symmetry condition, then we use biased formula.

          if bcl == -1 then
            y[1] = (1/(deltax^2))*(2*u[2] - 2*u[1]);
          else
            y[1] = (1/(deltax^2))*(u[3] - 2*u[2] + u[1]);
          end if;

          for i in 2:n-1 loop
            y[i] = (1/(deltax^2))*(u[i+1] - 2*u[i] + u[i-1]);
          end for;

          if bcr == -1 then
            y[n] = (1/(deltax^2))*(2*u[n-1] - 2*u[n]);
          else
            y[n] = (1/(deltax^2))*(u[n] - 2*u[n-1] + u[n-2]);
          end if;

          annotation (Icon(Text(
                extent=[-66,48; 66,-44],
                style(color=3, rgbcolor={0,0,255}),
                string="uxx_CD2B2")), Diagram,
            Documentation(info="<html>
<p>
Implements the second-order central difference approximation to the second-order space derivative.
</p>

<img align=middle src=\"..\\Images\\f16.png\">

<p>
For the boundary points second-order biased approximation is used.
</p>

<img align=middle src=\"..\\Images\\f14.png\">
<img align=middle src=\"..\\Images\\f15.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end u_xxCD2B2;

        block u_xxCD4B4
          extends PDE.Icons.BlockIcon;

        outer Integer n;
        outer Integer bcl;
        outer Integer bcr;

        parameter Real deltax = 1/(n-1);
        equation
          //bcl = -1 means: symmetry condition. If we do not have symmetry condition, then we use biased formula.

          if bcl == -1 then
            y[1] = (1/(12*deltax^2))*(-2*u[3] + 32*u[2] - 30*u[1]);
            y[2] = (1/(12*deltax^2))*(-u[4] + 16*u[3] - 31*u[2] + 16*u[1]);
          else
            y[1] = (1/(12*deltax^2))*(11*u[5] - 56*u[4] + 114*u[3] - 104*u[2] + 35*u[1]);
            y[2] = (1/(12*deltax^2))*(-u[5] + 4*u[4] + 6*u[3] - 20*u[2] + 11*u[1]);
          end if;

          for i in 3:n-2 loop
            y[i] = (1/(12*deltax^2))*(-u[i+2] + 16*u[i+1] - 30*u[i] + 16*u[i-1] - u[i-2]);
          end for;

          if bcr == -1 then
            y[n-1] = (1/(12*deltax^2))*(16*u[n] - 31*u[n-1] + 16*u[n-2] - u[n-3]);
            y[n] = (1/(12*deltax^2))*(-30*u[n] + 32*u[n-1] - 2*u[n-2]);
          else
            y[n-1] = (1/(12*deltax^2))*(11*u[n] - 20*u[n-1] + 6*u[n-2] + 4*u[n-3] - u[n-4]);
            y[n] = (1/(12*deltax^2))*(35*u[n] - 104*u[n-1] + 114*u[n-2] - 56*u[n-3] + 11*u[n-4]);
          end if;

        public
          Modelica.Blocks.Interfaces.RealInput u[n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-72,46; 70,-44],
                style(color=3, rgbcolor={0,0,255}),
                string="u_xxCD4B4")),
            Documentation(info="<html>
<p>
Implements the second-order central difference approximation to the second-order space derivative.
</p>

<img align=middle src=\"..\\Images\\f7.png\">

<p>
For the boundary points fourth-order biased approximation is used.
</p>

<img align=middle src=\"..\\Images\\f8.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end u_xxCD4B4;

        annotation (Documentation(info="<html>
<p>
This package contains blocks for the computation of spatial derivative schemes of different order.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end SDInterfaces;
    end SpaceDerivative;

    package Integrator
      block UniversalIntegrator
        extends PDE.Icons.BlockIcon1;

      outer PDE.World.worldModel worldModel1;
      parameter Integer n = worldModel1.n;

      parameter Integer vb = 1 "|Unknowns| The left most unknown";
      parameter Integer ve = worldModel1.n "|Unknowns| The right most unknown";

      parameter Integer icb = 1
          "|Initial Condition| Begin of the initial condition";
      parameter Integer ice = worldModel1.n
          "|Initial Condition| End of the initial condition";

      parameter Integer bcl = 0
          "|Boundary Conditions| Boundary condition at the left (0: no; 1: yes)";
      parameter Integer bcr = 0
          "|Boundary Conditions| Boundary condition at the right (0: no; 1: yes)";

      Real f[n];

      equation
        y = f;

        if bcl == 1 then
          f[1] = u2;
        end if;
        if bcr == 1 then
          f[n] = u3;
        end if;

        for i in vb:ve loop
          der(f[i]) = u[i];
        end for;

      initial equation

        for i in icb:ice loop
          f[i] = u1[i];
        end for;

        annotation (Icon(
            Text(
              extent=[-14,126; 16,80],
              style(color=3, rgbcolor={0,0,255}),
              string="%name"),
            Text(
              extent=[-98,-2; -68,-20],
              style(color=3, rgbcolor={0,0,255}),
              string="IC"),
            Text(
              extent=[-96,-32; -60,-48],
              style(color=3, rgbcolor={0,0,255}),
              string="BCL"),
            Text(
              extent=[-98,-62; -54,-78],
              style(color=3, rgbcolor={0,0,255}),
              string="BCR"),
            Text(
              extent=[-96,72; -68,48],
              style(color=3, rgbcolor={0,0,255}),
              string="R"),
            Text(
              extent=[54,74; 94,48],
              style(color=3, rgbcolor={0,0,255}),
              string="Var")), Diagram,
          Documentation(info="<html>
<p>
The integrator block accepts the equations of the form
</p>
<p>
<img align=middle src=\"..\\Images\\intdoc.png\">
</p>
<p>
where <i>f</i> is a function of <i>u</i>, <i>u_xx</i>, and so on. Once the equation is transformed in this form <br>
we can construct the right part <i>f</i> and pass it to the input <i>R</i> of the <b>integrator</b> block. <br>
The unknown varible <i>u</i> that we need for the construction of the right part of the equation is <br>
provided by the <i>Var</i> output of the integrator block.
The initial condition is passed to the <i>IC</i> input, whereas the left and right boundary conditions are <br>
passed to the <i>BCL</i> and <i>BCR</i> inputs respectively.
</p>
<p>
In the <i>PDE->MOL->Examples</i> Package many examples are implemented that show how to use the integrator block
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      public
        Modelica.Blocks.Interfaces.RealInput u[worldModel1.n]
          annotation (extent=[-118,50; -100,74]);
        Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
          annotation (extent=[100,50; 120,70]);
        Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n]
          annotation (extent=[-118,-22; -100,2]);
        Modelica.Blocks.Interfaces.RealInput u2
          annotation (extent=[-118,-52; -100,-28]);
        Modelica.Blocks.Interfaces.RealInput u3
          annotation (extent=[-118,-82; -100,-58]);
      equation

      end UniversalIntegrator;
      annotation (Documentation(info="<html>
<p>
This package contains integrator block that implements Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end Integrator;

    package Examples

      package Wave
        model WaveEquation
          inner PDE.World.worldModel worldModel1(qss=1,
            u_xx=1,
            n=10)                      annotation (extent=[-100,88; -60,100]);
          annotation (Diagram(
              Text(
                extent=[44,16; 76,10],
                style(color=3, rgbcolor={0,0,255}),
                string="(v.y, product.u2)"),
              Text(
                extent=[48,72; 70,66],
                style(color=3, rgbcolor={0,0,255}),
                string="(v.y, u.u)"),
              Text(
                extent=[-18,88; 18,82],
                style(color=3, rgbcolor={0,0,255}),
                string="(u.y, u_xx.u)")),
                               Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Wave equation</font></h3>
<p>
Implements the wave equation
</p>

<img align=middle src=\"..\\Images\\w1.png\">

<p>
where c is a constant value. The initial conditions are
</p>

<img align=middle src=\"..\\Images\\w3.png\">


<p>
and boundary conditions are
</p>

<img align=middle src=\"..\\Images\\w4.png\">

<p>
Because the integrator block cannot accept the equation in this form, we
transform the PDE above into two first-order PDEs:
</p>

<img align=middle src=\"..\\Images\\w2.png\">

<p>
The first equation is implemented in <b>u</b> block, the second in <b>v</b> block.
</p>
The analytical solution of this problem is implemented in <b>WaveAnalytic</b> block.
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Sources.RealExpression BCLu
            annotation (extent=[0,-4; 20,16]);
          Modelica.Blocks.Sources.RealExpression BCLv
            annotation (extent=[0,-66; 20,-46]);
          Modelica.Blocks.Sources.RealExpression ICv[worldModel1.n]
            annotation (extent=[0,-52; 20,-32]);
          WaveIC waveIC annotation (extent=[0,20; 20,40]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[-40,20; -20,40]);
          Modelica.Blocks.Sources.RealExpression Friction[worldModel1.n](y=0.0)
            annotation (extent=[-78,28; -66,44]);
          Modelica.Blocks.Math.Add add[worldModel1.n]
            annotation (extent=[0,42; 20,62]);
          WaveAnalytic waveAnalytic annotation (extent=[0,-92; 20,-72]);
          Modelica.Blocks.Math.Product product1[worldModel1.n]
            annotation (extent=[-40,52; -20,72]);
          Modelica.Blocks.Sources.RealExpression alpha2[worldModel1.n](y=1.0)
            annotation (extent=[-78,60; -66,76]);
          PDE.MOL.SpaceDerivative.Derivatives.u_xx u_xx(
                                    bcr=-1)
                                           annotation (extent=[-80,46; -64,60]);
          Integrator.UniversalIntegrator u(
            vb=2,
            icb=2,
            bcl=1) annotation (extent=[40,20; 80,60]);
          Integrator.UniversalIntegrator v(
            vb=2,
            icb=2,
            bcl=1) annotation (extent=[40,-60; 80,-20]);
        equation
          connect(product.y, add.u2) annotation (points=[-19,30; -16,30; -16,46;
                -2,46],
                     style(color=74, rgbcolor={0,0,127}));
          connect(Friction.y, product.u1) annotation (points=[-65.4,36; -42,36],
              style(color=74, rgbcolor={0,0,127}));
          connect(alpha2.y, product1.u1) annotation (points=[-65.4,68; -42,68],
              style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, add.u1) annotation (points=[-19,62; -10,62; -10,
                58; -2,58],                                               style(
                color=74, rgbcolor={0,0,127}));
          connect(u_xx.y, product1.u2) annotation (points=[-63.2,53; -52.6,53;
                -52.6,56; -42,56], style(color=74, rgbcolor={0,0,127}));
          connect(u.y, u_xx.u) annotation (points=[82,52; 88,52; 88,80; -90,80;
                -90,53; -81.6,53], style(color=74, rgbcolor={0,0,127}));
          connect(v.y, product.u2) annotation (points=[82,-28; 92,-28; 92,18;
                -46,18; -46,24; -42,24], style(color=74, rgbcolor={0,0,127}));
          connect(add.y, v.u) annotation (points=[21,52; 28,52; 28,-27.6; 38.2,
                -27.6], style(color=74, rgbcolor={0,0,127}));
          connect(v.y, u.u) annotation (points=[82,-28; 92,-28; 92,64; 32,64;
                32,52.4; 38.2,52.4], style(color=74, rgbcolor={0,0,127}));
          connect(waveIC.y, u.u1) annotation (points=[21,30; 30,30; 30,38; 38.2,
                38], style(color=74, rgbcolor={0,0,127}));
          connect(BCLu.y, u.u2) annotation (points=[21,6; 32,6; 32,32; 38.2,32],
              style(color=74, rgbcolor={0,0,127}));
          connect(ICv.y, v.u1) annotation (points=[21,-42; 38.2,-42], style(
                color=74, rgbcolor={0,0,127}));
          connect(BCLv.y, v.u2) annotation (points=[21,-56; 30,-56; 30,-48;
                38.2,-48], style(color=74, rgbcolor={0,0,127}));
        end WaveEquation;

        block WaveIC
          extends PDE.Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          inner parameter Integer n = worldModel1.n;
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          Modelica.Blocks.Math.Add add[worldModel1.n](k2=-1)
            annotation (extent=[-40,0; -20,20]);
          Modelica.Blocks.Sources.IntegerExpression integerExpression[worldModel1.n](y=1)
                    annotation (extent=[-80,-20; -60,0]);
          annotation (Diagram, Icon(Text(
                extent=[-58,44; 64,-44],
                style(color=3, rgbcolor={0,0,255}),
                string="WIC")),
            Documentation(info="<html>
<p>
Implements the initial condition for the first equation of the wave equation system
</p>

<img align=middle src=\"..\\Images\\w5.png\">


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[0,-10; 20,10]);
          Modelica.Blocks.Sources.RealExpression realExpression[worldModel1.n](y=3.14)
            annotation (extent=[-40,-40; -20,-20]);
          Modelica.Blocks.Sources.RealExpression realExpression1[worldModel1.n](y=2*(
                worldModel1.n - 1))
                               annotation (extent=[0,-40; 20,-20]);
          Modelica.Blocks.Math.Division division[worldModel1.n]
            annotation (extent=[50,-8; 66,8]);
          Modelica.Blocks.Sources.Constant index_i[worldModel1.n](k=1:worldModel1.
                n) annotation (extent=[-80,20; -60,40]);
          Modelica.Blocks.Math.Sin sin[worldModel1.n]
            annotation (extent=[76,-8; 92,8]);
        equation
          connect(realExpression.y, product.u2) annotation (points=[-19,-30; -10,
                -30; -10,-6; -2,-6],
                                  style(color=74, rgbcolor={0,0,127}));
          connect(add.y, product.u1) annotation (points=[-19,10; -10,10; -10,6;
                -2,6],
                     style(color=74, rgbcolor={0,0,127}));
          connect(product.y, division.u1) annotation (points=[21,0; 34,0; 34,
                4.8; 48.4,4.8],
                       style(color=74, rgbcolor={0,0,127}));
          connect(index_i.y, add.u1) annotation (points=[-59,30; -50,30; -50,16; -42,16],
              style(color=74, rgbcolor={0,0,127}));
          connect(realExpression1.y, division.u2) annotation (points=[21,-30;
                40,-30; 40,-4.8; 48.4,-4.8],
                                    style(color=74, rgbcolor={0,0,127}));
          connect(integerExpression.y, add.u2) annotation (points=[-59,-10; -50,
                -10; -50,4; -42,4], style(color=45, rgbcolor={255,127,0}));
          connect(division.y, sin.u) annotation (points=[66.8,0; 74.4,0], style(
                color=74, rgbcolor={0,0,127}));
          connect(sin.y, y) annotation (points=[92.8,0; 110,0], style(color=74,
                rgbcolor={0,0,127}));
        end WaveIC;

        block WaveAnalytic
          extends Icons.BlockIcon4;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
        equation
          for i in 1:n loop
            y[i] = 0.5*sin((3.14/2)*(((i-1)/(n-1)) - time)) + 0.5*sin((3.14/2)*(((i-1)/(n-1)) + time));
          end for;
          annotation (Diagram, Icon(Text(
                extent=[-54,28; 56,-24],
                style(color=3, rgbcolor={0,0,255}),
                string="WAN")),
            Documentation(info="<html>
<p>
Implements the analytical solution of the wave equation
</p>

<img align=middle src=\"..\\Images\\w6.png\">


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end WaveAnalytic;

        annotation (Documentation(info="<html>
<p>
This package contains wave equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Wave;

      package VibratingString
        model VibratingString
          inner PDE.World.worldModel worldModel1(qss=1, n=40)
                                       annotation (extent=[-100,88; -60,100]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Vibrating string equation</font></h3>
<p>
Implements the vibrating string equation. The finite vibrating string of length L (L = 1 in this problem) is
described by the wave equation together with special boundary and initial conditions. The complete problem is
</p>

<img align=middle src=\"..\\Images\\vs.png\">

<p>
where c is a constant value.
</p>


<p>
The analytical solution of this problem is implemented in <b>VibratingStringAnalytical</b> block.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Sources.RealExpression BCv
            annotation (extent=[0,-44; 20,-24]);
          Modelica.Blocks.Sources.RealExpression ICv[worldModel1.n]
            annotation (extent=[0,-28; 20,-8]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[-40,20; -20,40]);
          Modelica.Blocks.Sources.RealExpression Friction[worldModel1.n](y=0.0)
            annotation (extent=[-76,26; -62,46]);
          Modelica.Blocks.Math.Add add[worldModel1.n]
            annotation (extent=[0,40; 20,60]);
          Modelica.Blocks.Math.Product product1[worldModel1.n]
            annotation (extent=[-40,46; -20,66]);
          Modelica.Blocks.Sources.RealExpression alpha2[worldModel1.n](y=1.0)
            annotation (extent=[-58,54; -48,70]);
          PDE.MOL.Examples.VibratingString.VSIC vSIC
                    annotation (extent=[-78,0; -60,18]);
          PDE.MOL.Examples.VibratingString.VibratingStringAnalytical
            vibratingStringAnalytical
            annotation (extent=[0,-70; 20,-50]);
          Integrator.UniversalIntegrator u(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1) annotation (extent=[40,20; 80,60]);
          Integrator.UniversalIntegrator v(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1) annotation (extent=[40,-40; 80,0]);
          Modelica.Blocks.Sources.RealExpression BCu
            annotation (extent=[0,16; 20,36]);
          PDE.MOL.SpaceDerivative.Derivatives.u_xx u_xx
                                    annotation (extent=[-78,48; -60,60]);
        equation
          connect(product.y, add.u2) annotation (points=[-19,30; -16,30; -16,44; -2,
                44], style(color=74, rgbcolor={0,0,127}));
          connect(Friction.y, product.u1) annotation (points=[-61.3,36; -42,
                36],
              style(color=74, rgbcolor={0,0,127}));
          connect(alpha2.y, product1.u1) annotation (points=[-47.5,62; -42,62],
              style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, add.u1) annotation (points=[-19,56; -2,56],
              style(color=74, rgbcolor={0,0,127}));
          connect(BCv.y, v.u3) annotation (points=[21,-34; 38.2,-34], style(
                color=74, rgbcolor={0,0,127}));
          connect(BCv.y, v.u2) annotation (points=[21,-34; 28,-34; 28,-28;
                38.2,-28], style(color=74, rgbcolor={0,0,127}));
          connect(ICv.y, v.u1) annotation (points=[21,-18; 28,-18; 28,-22;
                38.2,-22], style(color=74, rgbcolor={0,0,127}));
          connect(BCu.y, u.u3) annotation (points=[21,26; 38.2,26], style(
                color=74, rgbcolor={0,0,127}));
          connect(BCu.y, u.u2) annotation (points=[21,26; 28,26; 28,32; 38.2,
                32], style(color=74, rgbcolor={0,0,127}));
          connect(vSIC.y, u.u1) annotation (points=[-59.1,9; -8,9; -8,38;
                38.2,38], style(color=74, rgbcolor={0,0,127}));
          connect(v.y, u.u) annotation (points=[82,-8; 96,-8; 96,66; 32,66;
                32,52.4; 38.2,52.4], style(color=74, rgbcolor={0,0,127}));
          connect(add.y, v.u) annotation (points=[21,50; 26,50; 26,-7.6; 38.2,
                -7.6], style(color=74, rgbcolor={0,0,127}));
          connect(v.y, product.u2) annotation (points=[82,-8; 88,-8; 88,16;
                -48,16; -48,24; -42,24], style(color=74, rgbcolor={0,0,127}));
          connect(u_xx.y, product1.u2) annotation (points=[-59.1,54; -52,54;
                -52,50; -42,50],
                         style(color=74, rgbcolor={0,0,127}));
          connect(u.y, u_xx.u) annotation (points=[82,52; 90,52; 90,76; -88,
                76; -88,54; -79.8,54], style(color=74, rgbcolor={0,0,127}));
        end VibratingString;

        block VSIC
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

        protected
                  Real pi = 3.14;

        equation
          for i in 1:n loop
            y[i] = sin(pi*(i-1)/(n-1)) + 0.5*sin(3*pi*(i-1)/(n-1));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Text(
                extent=[-42,26; 42,-26],
                style(color=3, rgbcolor={0,0,255}),
                string="VSIC")), Diagram,
            Documentation(info="<html>
<p>
Implements the initial condition
</p>

<img align=middle src=\"..\\Images\\vs2.png\">

<p>
for the <b>u</b> block of vibrating string equation.
</p>



</pre>
<p><b>Release Notes: </b></p>
</html>"));
        end VSIC;

        block VibratingStringAnalytical
          extends Icons.BlockIcon4;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Real alpha = 1.0;

        protected
                  Real pi = 3.14;

        equation
          for i in 1:n loop
            y[i] = (sin(pi*(i-1)/(n-1)))*(cos(pi*alpha*time)) + 0.5*(sin(3*pi*(i-1)/(n-1)))*(cos(3*pi*alpha*time));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Documentation(info="<html>
<p>
Implements the analytical solution of the vibrating string equation
</p>

<img align=middle src=\"..\\Images\\vs1.png\">


</pre>
<p><b>Release Notes: </b></p>
</html>"),
         Icon(Text(
                extent=[-50,40; 48,-38],
                style(color=3, rgbcolor={0,0,255}),
                string="VSA")));
        end VibratingStringAnalytical;
        annotation (Documentation(info="<html>
<p>
This package contains vibrating string equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end VibratingString;

      package Diffusion
        model DiffusionEquation

          PDE.MOL.Integrator.UniversalIntegrator Diffusion(
            vb=2,
            icb=2,
            bcl=1,
            n=worldModel1.n,
            ve=worldModel1.n - 1,
            ice=worldModel1.n - 1,
            bcr=1) annotation (extent=[20,0; 62,40]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Diffusion equation</font></h3>
<p>
Implements the diffusion equation
</p>

<img align=middle src=\"..\\Images\\d4.png\">


<p>
where sigma is a constant value. The initial condition is
</p>

<img align=middle src=\"..\\Images\\d1.png\">

<p>
and boundary conditions are
</p>

<img align=middle src=\"..\\Images\\d2.png\">

<p>
The analytical solution of this problem is implemented in <b>DiffusionAnalytic</b> block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          inner World.worldModel worldModel1(n=40)
            annotation (extent=[-100,88; -60,100]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[-18,34; -2,50]);
          Modelica.Blocks.Sources.RealExpression sigma[worldModel1.n](y=0.01)
                      annotation (extent=[-60,18; -40,38]);
          DIC1 dIC1_1 annotation (extent=[-18,10; -2,26]);
          DiffusionAnalytic diffusionAnalytic(alpha=0.1)
            annotation (extent=[-20,-40; 0,-20]);
          Modelica.Blocks.Sources.RealExpression BC
            annotation (extent=[-20,-16; 0,4]);
          PDE.MOL.SpaceDerivative.Derivatives.u_xx u_xx(
                                    bcr=-1) annotation (extent=[-60,40; -40,60]);
        equation
          connect(product.y, Diffusion.u)           annotation (points=[-1.2,42;
                8,42; 8,32.4; 18.11,32.4],   style(color=74, rgbcolor={0,0,127}));
          connect(dIC1_1.y, Diffusion.u1) annotation (points=[-1.2,18; 18.11,18], style(
                color=74, rgbcolor={0,0,127}));
          connect(BC.y, Diffusion.u2) annotation (points=[1,-6; 10,-6; 10,12;
                18.11,12], style(color=74, rgbcolor={0,0,127}));
          connect(BC.y, Diffusion.u3) annotation (points=[1,-6; 10,-6; 10,6;
                18.11,6], style(color=74, rgbcolor={0,0,127}));
          connect(sigma.y, product.u2) annotation (points=[-39,28; -28,28; -28,
                37.2; -19.6,37.2], style(color=74, rgbcolor={0,0,127}));
          connect(u_xx.y, product.u1) annotation (points=[-39,50; -28,50; -28,
                46.8; -19.6,46.8], style(color=74, rgbcolor={0,0,127}));
          connect(Diffusion.y, u_xx.u) annotation (points=[64.1,32; 72,32; 72,
                66; -70,66; -70,50; -62,50], style(color=74, rgbcolor={0,0,127}));
        end DiffusionEquation;

        block DiffusionIC
          extends PDE.Icons.BlockIcon;
          outer PDE.World.worldModel worldModel1;
          inner parameter Integer n = worldModel1.n;
          annotation (Diagram, Icon(Text(
                extent=[-50,42; 56,-40],
                style(color=3, rgbcolor={0,0,255}),
                string="DIC")),
            Documentation(info="<html>
<p>
Implements the initial condition cos(x) of the diffusion equation
</p>



</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Add add[worldModel1.n]
                                       annotation (extent=[-22,22; -2,42]);
          Modelica.Blocks.Sources.IntegerExpression integerExpression[worldModel1.n](
              y=-1)
            annotation (extent=[-80,6; -60,26]);
          Modelica.Blocks.Math.Division division[worldModel1.n]
                                                 annotation (extent=[20,0; 40,20]);
          Modelica.Blocks.Sources.IntegerExpression integerExpression1[worldModel1.n](y=n - 1)
            annotation (extent=[-20,-26; 0,-6]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
                                               annotation (extent=[60,-20; 80,0]);
          Modelica.Blocks.Sources.RealExpression realExpression[worldModel1.n](y=1.0)
            annotation (extent=[20,-36; 40,-16]);
          Modelica.Blocks.Math.Cos cos[worldModel1.n]
                                       annotation (extent=[60,-80; 80,-60]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          Modelica.Blocks.Sources.Constant index_i[worldModel1.n](k=1:worldModel1.
                n) annotation (extent=[-80,40; -60,60]);
        equation
          connect(integerExpression.y, add.u2) annotation (points=[-59,16; -42,16;
                -42,26; -24,26], style(color=45, rgbcolor={255,127,0}));
          connect(add.y, division.u1) annotation (points=[-1,32; 8,32; 8,16; 18,16],
              style(color=74, rgbcolor={0,0,127}));
          connect(integerExpression1.y, division.u2) annotation (points=[1,-16; 8,-16;
                8,4; 18,4], style(color=45, rgbcolor={255,127,0}));
          connect(division.y, product.u1) annotation (points=[41,10; 50,10; 50,-4; 58,
                -4], style(color=74, rgbcolor={0,0,127}));
          connect(realExpression.y, product.u2) annotation (points=[41,-26; 50,-26;
                50,-16; 58,-16], style(color=74, rgbcolor={0,0,127}));
          connect(cos.y, y) annotation (points=[81,-70; 92,-70; 92,0; 110,0], style(
                color=74, rgbcolor={0,0,127}));
          connect(index_i.y, add.u1) annotation (points=[-59,50; -42,50; -42,38;
                -24,38], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, cos.u) annotation (points=[81,-10; 86,-10; 86,-48;
                50,-48; 50,-70; 58,-70], style(color=74, rgbcolor={0,0,127}));
        end DiffusionIC;

        block DIC1
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

        equation
          for i in 1:n loop
            y[i] = sin(3.14*((i-1)/n)) + 0.5*sin(3*3.14*((i-1)/n));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Documentation(info="<html>
<p>
Implements the initial condition of the diffusion equation
</p>

<img align=middle src=\"..\\Images\\d1.png\">



</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"), Icon(Text(
                extent=[-60,42; 64,-40],
                style(color=3, rgbcolor={0,0,255}),
                string="DIC")));
        end DIC1;

        block DiffusionAnalytic
          extends Icons.BlockIcon4;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Real alpha = 0.0;

        equation
          for i in 1:n loop
            y[i] = exp(-((3.14*alpha)^2)*time)*sin(3.14*((i-1)/n)) + 0.5*exp(-((3*3.14*alpha)^2)*time)*sin(3*3.14*((i-1)/n));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Documentation(info="<html>
<p>
Implements the analytical solution of the diffusion equation
</p>

<img align=middle src=\"..\\Images\\d3.png\">



</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"), Icon(Text(
                extent=[-60,46; 58,-44],
                style(
                  color=3,
                  rgbcolor={0,0,255},
                  thickness=4,
                  fillColor=1,
                  rgbfillColor={255,0,0},
                  fillPattern=1),
                string="DAN")));
        end DiffusionAnalytic;
        annotation (Documentation(info="<html>
<p>
This package contains diffusion equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Diffusion;

      package Advection
        model AdvectionEquation

          inner World.worldModel worldModel1(qss=2, n=10)
            annotation (extent=[-100,88; -60,100]);
          PDE.MOL.Integrator.UniversalIntegrator Advection(
            ve=worldModel1.n,
            ice=worldModel1.n,
            bcr=0,
            vb=2,
            icb=2,
            bcl=1) annotation (extent=[20,0; 60,40]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Advection equation</font></h3>
<p>
Implements the linear advection equation
</p>

<img align=middle src=\"..\\Images\\a1.png\">

<p>
where c is a constant value. The initial condition is
</p>

<img align=middle src=\"..\\Images\\a3.png\">

<p>
and boundary condition at the left is
</p>

<img align=middle src=\"..\\Images\\a4.png\">

<p>
The analytical solution of this problem is implemented in <b>AdvectionAnalytic</b> block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[-18,32; -2,48]);
          Modelica.Blocks.Sources.RealExpression BCL(y=cos(-0.1*time))
            annotation (extent=[-20,-12; 0,8]);
          AdvectionAnalytic advectionAnalytic
            annotation (extent=[-20,-40; 0,-20]);
          Modelica.Blocks.Sources.RealExpression Speed[worldModel1.n](y=-0.1)
            annotation (extent=[-60,16; -40,36]);
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator
            annotation (extent=[-60,40; -40,60]);
          Diffusion.DiffusionIC diffusionIC annotation (extent=[-20,10; 0,26]);
        equation
          connect(BCL.y, Advection.u2)           annotation (points=[1,-2; 10,
                -2; 10,12; 18.2,12],  style(color=74, rgbcolor={0,0,127}));
          connect(diffusionIC.y, Advection.u1)           annotation (points=[1,18;
                18.2,18], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, Advection.u)           annotation (points=[-1.2,40;
                10,40; 10,32.4; 18.2,32.4], style(color=74, rgbcolor={0,0,127}));
          connect(Speed.y, product.u2) annotation (points=[-39,26; -28,26; -28,
                35.2; -19.6,35.2], style(color=74, rgbcolor={0,0,127}));
          connect(derivator.y, product.u1) annotation (points=[-38.9,50; -28,50;
                -28,44.8; -19.6,44.8], style(color=74, rgbcolor={0,0,127}));
          connect(Advection.y, derivator.u)           annotation (points=[62,32;
                68,32; 68,66; -70,66; -70,50; -62,50], style(color=74, rgbcolor=
                 {0,0,127}));
        end AdvectionEquation;

        block AdvectionAnalytic
          extends Icons.BlockIcon4;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Real speed = 0.1;

        equation
          for i in 1:n loop
            //y[i+1] = 1 - abs((i)/n);
            y[i] = cos(((i-1)/(n-1)) - speed*time);
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Documentation(info="<html>
<p>
Implements the analytical solution of the advection equation
</p>

<img align=middle src=\"..\\Images\\a2.png\">


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"),  Icon(Text(
                extent=[-62,34; 58,-30],
                style(
                  color=3,
                  rgbcolor={0,0,255},
                  thickness=4,
                  fillColor=1,
                  rgbfillColor={255,0,0},
                  fillPattern=1),
                string="AAN")));
        end AdvectionAnalytic;

        annotation (Documentation(info="<html>
<p>
This package contains advection equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Advection;

      package EulerSystem
        model ShockWaveEquation
          inner World.worldModel worldModel1(                          n=10, qss=1)
            annotation (extent=[-100,88; -60,100]);
          PDE.MOL.Integrator.UniversalIntegrator Density(
            ve=worldModel1.n - 1,
            ice=worldModel1.n - 1,
            bcr=1)                 annotation (extent=[-42,42; -10,66]);
          PDE.MOL.Integrator.UniversalIntegrator Velocity(
            vb=2,
            icb=2,
            bcl=1) annotation (extent=[-42,-8; -8,14]);
          PDE.MOL.Integrator.UniversalIntegrator Pressure(
            ve=worldModel1.n - 1,
            ice=worldModel1.n - 1,
            bcr=1)                 annotation (extent=[-42,-72; -8,-50]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[26,38; 38,50]);
          Modelica.Blocks.Math.Product product1[worldModel1.n]
            annotation (extent=[26,58; 38,70]);
          Modelica.Blocks.Math.Gain gain[worldModel1.n](k=-1)
            annotation (extent=[46,60; 52,68]);
          Modelica.Blocks.Math.Gain gain1[worldModel1.n](k=-1)
            annotation (extent=[46,40; 52,48]);
          Modelica.Blocks.Math.Add add[worldModel1.n]
            annotation (extent=[66,48; 78,60]);
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator(
                                                      bcl=0, bcr=-1)
                                              annotation (extent=[6,76; 16,86]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Euler equations</font></h3>
<p>
Implements the Euler system of equations
</p>
<p>
<img align=middle src=\"..\\Images\\sw1.png\">
</p>

<p>
where
</p>

<img align=middle src=\"..\\Images\\sw3.png\">

<p>
The initial conditions are
</p>

<img align=middle src=\"..\\Images\\sw4.png\">
<p>
where <i>rho<sub>i</sub></i> is determined by the equation of state for ideal gases:
</p>

<img align=middle src=\"..\\Images\\sw5.png\">

<p>
where T = 300K, R = 8.314 J/K*mole and M<sub>air</sub> = 28.96 g/mole
</p>


<p>
The boundary conditions are
</p>

<img align=middle src=\"..\\Images\\sw6.png\">




</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator1(
                                                       bcl=-1)
                                               annotation (extent=[6,16; 16,26]);
          Modelica.Blocks.Math.Division division[worldModel1.n]
            annotation (extent=[-22,-30; -12,-20]);
          Modelica.Blocks.Sources.RealExpression realExpr[worldModel1.n](y=1)
            annotation (extent=[-40,-30; -32,-12]);
          Modelica.Blocks.Math.Product product2[worldModel1.n]
            annotation (extent=[26,-6; 38,6]);
          Modelica.Blocks.Math.Product product3[worldModel1.n]
            annotation (extent=[6,-38; 18,-26]);
          Modelica.Blocks.Math.Add add1[worldModel1.n]
            annotation (extent=[66,-22; 78,-10]);
          Modelica.Blocks.Math.Gain gain2[worldModel1.n](k=-1)
            annotation (extent=[46,-4; 52,4]);
          Modelica.Blocks.Math.Gain gain3[worldModel1.n](k=-1)
            annotation (extent=[46,-36; 52,-28]);
          Modelica.Blocks.Math.Product product4[worldModel1.n]
            annotation (extent=[6,-60; 18,-48]);
          Modelica.Blocks.Math.Product product5[worldModel1.n]
            annotation (extent=[-2,-90; 10,-78]);
          Modelica.Blocks.Math.Add add2[worldModel1.n]
            annotation (extent=[66,-74; 78,-62]);
          Modelica.Blocks.Math.Gain gain4[worldModel1.n](k=-1)
            annotation (extent=[46,-58; 52,-50]);
          Modelica.Blocks.Math.Gain gain5[worldModel1.n](k=-1)
            annotation (extent=[46,-86; 52,-78]);
          Modelica.Blocks.Math.Product product6[worldModel1.n]
            annotation (extent=[26,-88; 38,-76]);
          Modelica.Blocks.Sources.RealExpression gamma[worldModel1.n](y=1.4)
            annotation (extent=[-26,-96; -18,-78]);
          Modelica.Blocks.Sources.RealExpression ICdensity[worldModel1.n](y=1.0)
            annotation (extent=[-84,50; -76,66]);
          Modelica.Blocks.Sources.RealExpression ICvelocity[worldModel1.n](y=0.0)
            annotation (extent=[-84,-10; -76,8]);
          Modelica.Blocks.Sources.RealExpression ICpressure[worldModel1.n](y=1.0)
                        annotation (extent=[-84,-72; -74,-54]);
          Modelica.Blocks.Sources.RealExpression BCR(y=0.125)
            annotation (extent=[-84,36; -76,52]);
          Modelica.Blocks.Sources.RealExpression BCL(y=0.0)
            annotation (extent=[-84,-24; -76,-6]);
          Modelica.Blocks.Sources.RealExpression BCRp(y=0.1)
            annotation (extent=[-84,-86; -74,-68]);
          a a1 annotation (extent=[-20,-44; -10,-34]);
        equation
          connect(Density.y, derivator.u) annotation (points=[-8.4,61.2; -2,
                61.2; -2,81; 5,81],
                           style(color=74, rgbcolor={0,0,127}));
          connect(derivator.y, product1.u1) annotation (points=[16.55,81; 20,81;
                20,67.6; 24.8,67.6], style(color=74, rgbcolor={0,0,127}));
          connect(Velocity.y, derivator1.u) annotation (points=[-6.3,9.6; 2,9.6;
                2,21; 5,21],
                       style(color=74, rgbcolor={0,0,127}));
          connect(Density.y, product.u1) annotation (points=[-8.4,61.2; -2,61.2;
                -2,47.6; 24.8,47.6],
                                  style(color=74, rgbcolor={0,0,127}));
          connect(Velocity.y, product1.u2) annotation (points=[-6.3,9.6; -4,9.6;
                -4,38; 6,38; 6,60.4; 24.8,60.4],
                                              style(color=74, rgbcolor={0,0,127}));
          connect(derivator1.y, product.u2) annotation (points=[16.55,21; 20,21;
                20,40.4; 24.8,40.4], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, gain.u) annotation (points=[38.6,64; 45.4,64],
              style(color=74, rgbcolor={0,0,127}));
          connect(product.y, gain1.u) annotation (points=[38.6,44; 45.4,44],
              style(color=74, rgbcolor={0,0,127}));
          connect(gain.y, add.u1) annotation (points=[52.3,64; 58,64; 58,57.6;
                64.8,57.6], style(color=74, rgbcolor={0,0,127}));
          connect(gain1.y, add.u2) annotation (points=[52.3,44; 58,44; 58,50.4;
                64.8,50.4], style(color=74, rgbcolor={0,0,127}));
          connect(add.y, Density.u) annotation (points=[78.6,54; 86,54; 86,90;
                -54,90; -54,61.44; -43.44,61.44],
                                                style(color=74, rgbcolor={0,0,127}));
          connect(realExpr.y, division.u1) annotation (points=[-31.6,-21; -27.8,
                -21; -27.8,-22; -23,-22], style(color=74, rgbcolor={0,0,127}));
          connect(division.y, product3.u1) annotation (points=[-11.5,-25; -3.75,
                -25; -3.75,-28.4; 4.8,-28.4], style(color=74, rgbcolor={0,0,127}));
          connect(Velocity.y, product2.u2) annotation (points=[-6.3,9.6; -4,9.6;
                -4,-3.6; 24.8,-3.6],
                                  style(color=74, rgbcolor={0,0,127}));
          connect(derivator1.y, product2.u1) annotation (points=[16.55,21; 20,
                21; 20,3.6; 24.8,3.6],
                                   style(color=74, rgbcolor={0,0,127}));
          connect(product3.y, gain3.u) annotation (points=[18.6,-32; 45.4,-32],
              style(color=74, rgbcolor={0,0,127}));
          connect(product2.y, gain2.u) annotation (points=[38.6,0; 45.4,0], style(
                color=74, rgbcolor={0,0,127}));
          connect(gain2.y, add1.u1) annotation (points=[52.3,0; 58,0; 58,-12.4;
                64.8,-12.4], style(color=74, rgbcolor={0,0,127}));
          connect(gain3.y, add1.u2) annotation (points=[52.3,-32; 58,-32; 58,
                -19.6; 64.8,-19.6], style(color=74, rgbcolor={0,0,127}));
          connect(product4.y, gain4.u) annotation (points=[18.6,-54; 45.4,-54],
              style(color=74, rgbcolor={0,0,127}));
          connect(derivator1.y, product6.u1) annotation (points=[16.55,21; 20,21;
                20,-78.4; 24.8,-78.4], style(color=74, rgbcolor={0,0,127}));
          connect(product5.y, product6.u2) annotation (points=[10.6,-84; 17.7,-84;
                17.7,-85.6; 24.8,-85.6], style(color=74, rgbcolor={0,0,127}));
          connect(gamma.y, product5.u2) annotation (points=[-17.6,-87; -10.8,
                -87; -10.8,-87.6; -3.2,-87.6],
                                          style(color=74, rgbcolor={0,0,127}));
          connect(Pressure.y, product5.u1) annotation (points=[-6.3,-54.4; -6,
                -54.4; -6,-80.4; -3.2,-80.4],
                                       style(color=74, rgbcolor={0,0,127}));
          connect(product6.y, gain5.u) annotation (points=[38.6,-82; 45.4,-82],
              style(color=74, rgbcolor={0,0,127}));
          connect(gain4.y, add2.u1) annotation (points=[52.3,-54; 58,-54; 58,
                -64.4; 64.8,-64.4], style(color=74, rgbcolor={0,0,127}));
          connect(gain5.y, add2.u2) annotation (points=[52.3,-82; 58,-82; 58,
                -71.6; 64.8,-71.6], style(color=74, rgbcolor={0,0,127}));
          connect(BCR.y, Density.u3) annotation (points=[-75.6,44; -48,44; -48,
                45.6; -43.44,45.6],
                            style(color=74, rgbcolor={0,0,127}));
          connect(BCL.y, Velocity.u2) annotation (points=[-75.6,-15; -48,-15;
                -48,-1.4; -43.53,-1.4],
                                style(color=74, rgbcolor={0,0,127}));
          connect(BCRp.y, Pressure.u3) annotation (points=[-73.5,-77; -48,-77;
                -48,-68.7; -43.53,-68.7],
                                      style(color=74, rgbcolor={0,0,127}));
          connect(Density.y, division.u2) annotation (points=[-8.4,61.2; -6,
                61.2; -6,34; -92,34; -92,-28; -23,-28],
                                               style(color=74, rgbcolor={0,0,127}));
          connect(Velocity.y, product4.u2) annotation (points=[-6.3,9.6; -2,9.6;
                -2,-57.6; 4.8,-57.6],
                                   style(color=74, rgbcolor={0,0,127}));
          connect(add2.y, Pressure.u) annotation (points=[78.6,-68; 86,-68; 86,
                -96; -96,-96; -96,-54.18; -43.53,-54.18],
                                                        style(color=74, rgbcolor=
                  {0,0,127}));
          connect(ICpressure.y, Pressure.u1) annotation (points=[-73.5,-63;
                -57.75,-63; -57.75,-62.1; -43.53,-62.1],
                                                     style(color=74, rgbcolor={0,
                  0,127}));
          connect(add1.y, Velocity.u) annotation (points=[78.6,-16; 86,-16; 86,
                28; -52,28; -52,9.82; -43.53,9.82],
                                              style(color=74, rgbcolor={0,0,127}));
          connect(a1.y, product3.u2) annotation (points=[-9.5,-39; -6,-39; -6,
                -35.6; 4.8,-35.6],
                        style(color=74, rgbcolor={0,0,127}));
          connect(a1.y, product4.u1) annotation (points=[-9.5,-39; -6,-39; -6,
                -50.4; 4.8,-50.4], style(color=74, rgbcolor={0,0,127}));
          connect(Density.y, a1.u) annotation (points=[-8.4,61.2; -6,61.2; -6,
                34; -92,34; -92,-36; -21,-36],
                                           style(color=74, rgbcolor={0,0,127}));
          connect(Velocity.y, a1.u1) annotation (points=[-6.3,9.6; -4,9.6; -4,
                20; -94,20; -94,-39; -21,-39],
                                           style(color=74, rgbcolor={0,0,127}));
          connect(Pressure.y, a1.u2) annotation (points=[-6.3,-54.4; -4,-54.4;
                -4,-46; -24,-46; -24,-42; -21,-42],
                                                 style(color=74, rgbcolor={0,0,
                  127}));
          connect(ICdensity.y, Density.u1) annotation (points=[-75.6,58; -58,58;
                -58,52.8; -43.44,52.8],
                                    style(color=74, rgbcolor={0,0,127}));
          connect(ICvelocity.y, Velocity.u1) annotation (points=[-75.6,-1; -60,
                -1; -60,1.9; -43.53,1.9], style(color=74, rgbcolor={0,0,127}));
        end ShockWaveEquation;

        block Friction
          extends PDE.Icons.BlockIcon;

        outer Integer n;
        parameter Real alpha = 0.1;

          Modelica.Blocks.Interfaces.RealInput u[n]
            annotation (extent=[-140,20; -100,60]);
          Modelica.Blocks.Interfaces.RealInput u1[n]
            annotation (extent=[-140,-60; -100,-20]);
          Modelica.Blocks.Interfaces.RealOutput y[n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-44,34; 44,-32],
                style(color=3, rgbcolor={0,0,255}),
                string="Friction")),
            Documentation(info="<html>
<p>
Implements the frictional resistance of the Euler system of equations
</p>
<p>
<img align=middle src=\"..\\Images\\fb.png\">
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Product product[n]
            annotation (extent=[-60,-10; -40,10]);
          Modelica.Blocks.Math.Abs abs1[n]
            annotation (extent=[-60,-50; -40,-30]);
          Modelica.Blocks.Math.Product product1[n]
            annotation (extent=[-20,-30; 0,-10]);
          Modelica.Blocks.Math.Product product2[n]
            annotation (extent=[20,0; 40,20]);
          Modelica.Blocks.Sources.RealExpression param_alpha[n](y=0.1)
            annotation (extent=[-20,20; 0,40]);
          Modelica.Blocks.Math.Division division[n]
            annotation (extent=[60,-10; 80,10]);
          Modelica.Blocks.Sources.RealExpression delta_x[n](y=1/(
                n - 1))             annotation (extent=[20,-40; 40,-20]);
        equation
          connect(u, product.u1) annotation (points=[-120,40; -92,40; -92,6; -62,6],
              style(color=74, rgbcolor={0,0,127}));
          connect(u1, product.u2) annotation (points=[-120,-40; -92,-40; -92,-6; -62,
                -6], style(color=74, rgbcolor={0,0,127}));
          connect(u1, abs1.u) annotation (points=[-120,-40; -62,-40], style(color=74,
                rgbcolor={0,0,127}));
          connect(product.y, product1.u1) annotation (points=[-39,0; -32,0; -32,-14;
                -22,-14], style(color=74, rgbcolor={0,0,127}));
          connect(abs1.y, product1.u2) annotation (points=[-39,-40; -32,-40; -32,-26;
                -22,-26], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, product2.u2) annotation (points=[1,-20; 10,-20; 10,4;
                18,4], style(color=74, rgbcolor={0,0,127}));
          connect(param_alpha.y, product2.u1) annotation (points=[1,30; 10,30; 10,16;
                18,16], style(color=74, rgbcolor={0,0,127}));
          connect(product2.y, division.u1) annotation (points=[41,10; 50,10; 50,6; 58,
                6], style(color=74, rgbcolor={0,0,127}));
          connect(delta_x.y, division.u2) annotation (points=[41,-30; 50,-30; 50,-6;
                58,-6], style(color=74, rgbcolor={0,0,127}));
          connect(division.y, y)
            annotation (points=[81,0; 110,0], style(color=74, rgbcolor={0,0,127}));
        end Friction;

        block Viscosity
          extends PDE.Icons.BlockIcon;

          outer Integer n;

          Modelica.Blocks.Interfaces.RealInput u[n]
            annotation (extent=[-140,20; -100,60]);
          Modelica.Blocks.Interfaces.RealInput u1[n]
            annotation (extent=[-140,-60; -100,-20]);
          Modelica.Blocks.Interfaces.RealOutput y[n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-54,24; 52,-20],
                style(color=3, rgbcolor={0,0,255}),
                string="Viscosity")),
            Documentation(info="<html>
<p>
Implements the pseudo viscous pressure of the Euler system of equations
</p>
<p>
<img align=middle src=\"..\\Images\\vb.png\">
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Logical.Less less[n]
            annotation (extent=[-60,-50; -40,-30]);
          Modelica.Blocks.Sources.RealExpression realExpression[n]
            annotation (extent=[-84,-66; -70,-46]);
          Modelica.Blocks.Logical.Switch switch1[n]
            annotation (extent=[60,-10; 80,10]);
          Modelica.Blocks.Sources.RealExpression realExpression1[n]
            annotation (extent=[6,-34; 20,-14]);
          Modelica.Blocks.Sources.RealExpression param_beta[n](y=0.1)
            annotation (extent=[-62,76; -48,96]);
          Modelica.Blocks.Sources.RealExpression delta_x[n](y=1/(
                n - 1))             annotation (extent=[-80,60; -66,80]);
          Modelica.Blocks.Math.Product product[n]
            annotation (extent=[-32,24; -12,44]);
          Modelica.Blocks.Math.Product product1[n]
            annotation (extent=[-70,-12; -50,8]);
          Modelica.Blocks.Math.Product product2[n]
            annotation (extent=[-44,54; -24,74]);
          Modelica.Blocks.Math.Product product3[n]
            annotation (extent=[-8,70; 12,90]);
          Modelica.Blocks.Math.Product product4[n]
            annotation (extent=[24,40; 44,60]);
        equation
          connect(u1, less.u1) annotation (points=[-120,-40; -62,-40], style(color=74,
                rgbcolor={0,0,127}));
          connect(realExpression.y, less.u2) annotation (points=[-69.3,-56; -66,
                -56; -66,-48; -62,-48],
                                   style(color=74, rgbcolor={0,0,127}));
          connect(less.y, switch1.u2) annotation (points=[-39,-40; -30,-40; -30,0; 58,
                0], style(color=5, rgbcolor={255,0,255}));
          connect(realExpression1.y, switch1.u3) annotation (points=[20.7,-24;
                30,-24; 30,-8; 58,-8],
                               style(color=74, rgbcolor={0,0,127}));
          connect(u, product.u1) annotation (points=[-120,40; -34,40],
              style(color=74, rgbcolor={0,0,127}));
          connect(u1, product1.u1) annotation (points=[-120,-40; -98,-40; -98,4; -72,
                4], style(color=74, rgbcolor={0,0,127}));
          connect(u1, product1.u2) annotation (points=[-120,-40; -96,-40; -96,-8; -72,
                -8], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, product.u2) annotation (points=[-49,-2; -46,-2; -46,
                28; -34,28],
                         style(color=74, rgbcolor={0,0,127}));
          connect(delta_x.y, product2.u1) annotation (points=[-65.3,70; -46,70],
              style(color=74, rgbcolor={0,0,127}));
          connect(delta_x.y, product2.u2) annotation (points=[-65.3,70; -56,70;
                -56,58; -46,58],
                             style(color=74, rgbcolor={0,0,127}));
          connect(param_beta.y, product3.u1) annotation (points=[-47.3,86; -10,86],
              style(color=74, rgbcolor={0,0,127}));
          connect(product2.y, product3.u2) annotation (points=[-23,64; -18,64; -18,74;
                -10,74], style(color=74, rgbcolor={0,0,127}));
          connect(product3.y, product4.u1) annotation (points=[13,80; 18,80; 18,56;
                22,56], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, product4.u2) annotation (points=[-11,34; 2,34; 2,44;
                22,44],
                     style(color=74, rgbcolor={0,0,127}));
          connect(product4.y, switch1.u1) annotation (points=[45,50; 52,50; 52,8; 58,
                8], style(color=74, rgbcolor={0,0,127}));
          connect(switch1.y, y)
            annotation (points=[81,0; 110,0], style(color=74, rgbcolor={0,0,127}));
        end Viscosity;

        block a
          extends PDE.Icons.BlockIcon;
          outer PDE.World.worldModel worldModel1;
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n]
            annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n]
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Text(
                extent=[-52,44; 56,-32],
                style(color=3, rgbcolor={0,0,255}),
                string="a")), Diagram,
            Documentation(info="<html>
<p>
Implements
</p>
<p>
<img align=middle src=\"..\\Images\\ab.png\">
</p>
<p>
in the Euler system of equations
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator1
                               annotation (extent=[0,-10; 20,10]);
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator2
                               annotation (extent=[-44,-70; -24,-50]);
          Modelica.Blocks.Math.Add3 add3_1[worldModel1.n]
            annotation (extent=[60,-10; 80,10]);
          parameter Integer method = worldModel1.qss;
          inner parameter Integer n = worldModel1.n;
          Friction friction annotation (extent=[-40,60; -20,80]);
          Viscosity viscosity annotation (extent=[-40,14; -20,34]);
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator3
            annotation (extent=[-74,-30; -54,-10]);
        equation
          connect(u2, derivator2.u) annotation (points=[-120,-60; -46,-60],
                          style(color=74, rgbcolor={0,0,127}));
          connect(derivator1.y, add3_1.u2) annotation (points=[21.1,0; 58,0],
                       style(color=74, rgbcolor={0,0,127}));
          connect(add3_1.y, y)
            annotation (points=[81,0; 110,0], style(color=74, rgbcolor={0,0,127}));
          connect(derivator2.y, add3_1.u3) annotation (points=[-22.9,-60; 38,
                -60; 38,-8; 58,-8],
                        style(color=74, rgbcolor={0,0,127}));
          connect(u, friction.u) annotation (points=[-120,60; -88,60; -88,74;
                -42,74], style(color=74, rgbcolor={0,0,127}));
          connect(u1, friction.u1) annotation (points=[-120,0; -82,0; -82,66;
                -42,66], style(color=74, rgbcolor={0,0,127}));
          connect(viscosity.y, derivator1.u) annotation (points=[-19,24; -10,24;
                -10,0; -2,0],
              style(color=74, rgbcolor={0,0,127}));
          connect(u, viscosity.u) annotation (points=[-120,60; -88,60; -88,28;
                -42,28],style(color=74, rgbcolor={0,0,127}));
          connect(u1, derivator3.u) annotation (points=[-120,0; -82,0; -82,-20;
                -76,-20], style(color=74, rgbcolor={0,0,127}));
          connect(derivator3.y, viscosity.u1) annotation (points=[-52.9,-20;
                -48,-20; -48,20; -42,20], style(color=74, rgbcolor={0,0,127}));
          connect(friction.y, add3_1.u1) annotation (points=[-19,70; 38,70; 38,8; 58,
                8], style(color=74, rgbcolor={0,0,127}));
        end a;

        annotation (Documentation(info="<html>
<p>
This package contains Euler equations solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end EulerSystem;

      package Burger
        model BurgerEquation
          inner World.worldModel worldModel1(             qss=2, n=40)
            annotation (extent=[-20,88; 20,100]);
          Integrator.UniversalIntegrator Burger(
            vb=2,
            icb=2,
            bcl=1,
            ve=worldModel1.n,
            ice=worldModel1.n,
            bcr=0)
            annotation (extent=[-60,20; -20,60]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[0,46; 12,58]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Burger큦 equation</font></h3>
<p>
Implements the inviscid Burger큦 equation
</p>

<img align=middle src=\"..\\Images\\b1.png\">

<p>
The initial condition is
</p>

<img align=middle src=\"..\\Images\\b2.png\">

<p>
and boundary conditions are
</p>

<p>
<img align=middle src=\"..\\Images\\b3.png\">
</p>


<p>
The analytical solution of this problem is implemented in <b>BAN</b> block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Product product1[worldModel1.n]
            annotation (extent=[48,44; 58,54]);
          Modelica.Blocks.Sources.RealExpression const[worldModel1.n](y=-0.5)
            annotation (extent=[26,26; 34,42]);
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator
                                              annotation (extent=[24,46; 36,58]);
          BICx bICx annotation (extent=[-86,32; -74,44]);
          BAN bAN annotation (extent=[0,-20; 20,0]);
          Modelica.Blocks.Sources.RealExpression BC
            annotation (extent=[-90,10; -70,30]);
        equation
          connect(Burger.y, product.u1) annotation (points=[-18,52; -10,52; -10,
                55.6; -1.2,55.6], style(color=74, rgbcolor={0,0,127}));
          connect(Burger.y, product.u2) annotation (points=[-18,52; -10,52; -10,
                48.4; -1.2,48.4], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, derivator.u) annotation (points=[12.6,52; 22.8,52],
              style(color=74, rgbcolor={0,0,127}));
          connect(derivator.y, product1.u1) annotation (points=[36.66,52; 47,52],
              style(color=74, rgbcolor={0,0,127}));
          connect(const.y, product1.u2) annotation (points=[34.4,34; 42,34; 42,
                46; 47,46], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, Burger.u) annotation (points=[58.5,49; 66,49; 66,
                66; -70,66; -70,52.4; -61.8,52.4], style(color=74, rgbcolor={0,
                  0,127}));
          connect(bICx.y, Burger.u1) annotation (points=[-73.4,38; -61.8,38],
                               style(color=74, rgbcolor={0,0,127}));
          connect(BC.y, Burger.u2)             annotation (points=[-69,20; -66,
                20; -66,32; -61.8,32],
                                   style(color=74, rgbcolor={0,0,127}));
          connect(BC.y, Burger.u3)             annotation (points=[-69,20; -66,
                20; -66,26; -61.8,26],
                                   style(color=74, rgbcolor={0,0,127}));
          connect(Burger.y, bAN.u) annotation (points=[-18,52; -10,52; -10,-10;
                -2,-10], style(color=74, rgbcolor={0,0,127}));
        end BurgerEquation;

        block BICx
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

        equation
          for i in 1:n loop
            y[i] = (i-1)/(n-1);
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Documentation(info="<html>
<p>
Implements the initial condition of the inviscid Burger큦 equation
</p>

<img align=middle src=\"..\\Images\\b2.png\">


</pre>
<p><b>Release Notes: </b></p>
</html>"),  Icon(Text(
                extent=[-52,38; 48,-34],
                style(color=3, rgbcolor={0,0,255}),
                string="BIC")));
        end BICx;

        block BAN
          extends Icons.BlockIcon4;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n]
            annotation (extent=[-140,-20; -100,20]);
        equation
          for i in 1:n loop
            //y[i] = (i-1)/(n-1) - u[i]*time;
            y[i] = ((i-1)/(n-1))/(1 + time);
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Documentation(info="<html>
<p>
Implements the analytical solution of the inviscid Burger큦 equation
</p>

<img align=middle src=\"..\\Images\\b5.png\">
</pre>
<p><b>Release Notes: </b></p>
</html>"),  Icon(Text(
                extent=[-54,40; 60,-38],
                style(
                  color=3,
                  rgbcolor={0,0,255},
                  thickness=4,
                  fillColor=1,
                  rgbfillColor={255,0,0},
                  fillPattern=1),
                string="BAN")));
        end BAN;
        annotation (Documentation(info="<html>
<p>
This package contains Burger큦 equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Burger;

      package BuckleyLeverett
        model BuckleyLeverettEquation
          inner World.worldModel worldModel1(qss=2, n=10)
            annotation (extent=[-20,88; 20,100]);
          Integrator.UniversalIntegrator BuckleyLeverett(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1)
            annotation (extent=[-64,12; -30,44]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[-18,32; -8,42]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Buckley-Leverett equation</font></h3>
<p>
Implements the Buckley-Leverett equation
</p>

<img align=middle src=\"..\\Images\\bl1.png\">

<p>
The initial condition is
</p>

<img align=middle src=\"..\\Images\\bl2.png\">

<p>
and boundary conditions are
</p>

<p>
<img align=middle src=\"..\\Images\\bl3.png\">
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Product product1[worldModel1.n]
            annotation (extent=[0,28; 10,38]);
          Modelica.Blocks.Sources.RealExpression const[worldModel1.n](y=4.0)
            annotation (extent=[-18,16; -8,32]);
          Modelica.Blocks.Math.Division division[worldModel1.n]
            annotation (extent=[56,24; 66,34]);
          Modelica.Blocks.Math.Add add[worldModel1.n](k2=-1)
            annotation (extent=[0,0; 10,10]);
          Modelica.Blocks.Sources.RealExpression one[worldModel1.n](y=1.0)
            annotation (extent=[-18,0; -8,16]);
          Modelica.Blocks.Math.Product product2[worldModel1.n]
            annotation (extent=[20,0; 30,10]);
          Modelica.Blocks.Math.Add add1[worldModel1.n]
            annotation (extent=[38,14; 48,24]);
          PDE.MOL.SpaceDerivative.Derivatives.u_x derivator
                                              annotation (extent=[70,24; 80,34]);
          Modelica.Blocks.Math.Gain gain[worldModel1.n](k=-1)
            annotation (extent=[84,26; 88,32]);
          Modelica.Blocks.Sources.RealExpression IC[worldModel1.n](y=1.0)
            annotation (extent=[-88,26; -80,40]);
          Modelica.Blocks.Sources.RealExpression BCL
            annotation (extent=[-88,14; -80,28]);
          Modelica.Blocks.Sources.RealExpression BCR
            annotation (extent=[-88,2; -80,16]);
        equation
          connect(BuckleyLeverett.y, product.u1) annotation (points=[-28.3,37.6;
                -22,37.6; -22,40; -19,40],style(color=74, rgbcolor={0,0,127}));
          connect(BuckleyLeverett.y, product.u2) annotation (points=[-28.3,37.6;
                -22,37.6; -22,34; -19,34],style(color=74, rgbcolor={0,0,127}));
          connect(product.y, product1.u1) annotation (points=[-7.5,37; -3.75,37;
                -3.75,36; -1,36],
                                style(color=74, rgbcolor={0,0,127}));
          connect(const.y, product1.u2) annotation (points=[-7.5,24; -4,24; -4,
                30; -1,30],
                       style(color=74, rgbcolor={0,0,127}));
          connect(one.y, add.u1) annotation (points=[-7.5,8; -1,8],
                                                                  style(color=
                  74, rgbcolor={0,0,127}));
          connect(BuckleyLeverett.y, add.u2) annotation (points=[-28.3,37.6;
                -24,37.6; -24,2; -1,2],style(color=74, rgbcolor={0,0,127}));
          connect(add.y, product2.u1) annotation (points=[10.5,5; 14,5; 14,8;
                19,8], style(color=74, rgbcolor={0,0,127}));
          connect(add.y, product2.u2) annotation (points=[10.5,5; 14,5; 14,2;
                19,2], style(color=74, rgbcolor={0,0,127}));
          connect(product2.y, add1.u2) annotation (points=[30.5,5; 34,5; 34,16;
                37,16], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, add1.u1) annotation (points=[10.5,33; 34,33; 34,
                22; 37,22], style(color=74, rgbcolor={0,0,127}));
          connect(add1.y, division.u2) annotation (points=[48.5,19; 52,19; 52,
                26; 55,26], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, division.u1) annotation (points=[10.5,33; 50,33;
                50,32; 55,32], style(color=74, rgbcolor={0,0,127}));
          connect(division.y, derivator.u) annotation (points=[66.5,29; 69,29],
              style(color=74, rgbcolor={0,0,127}));
          connect(derivator.y, gain.u) annotation (points=[80.55,29; 83.6,29],
                                     style(color=74, rgbcolor={0,0,127}));
          connect(gain.y, BuckleyLeverett.u) annotation (points=[88.2,29; 90,29;
                90,50; -72,50; -72,37.92; -65.53,37.92],      style(color=74,
                rgbcolor={0,0,127}));
          connect(BCR.y, BuckleyLeverett.u3) annotation (points=[-79.6,9; -72,9;
                -72,16.8; -65.53,16.8], style(color=74, rgbcolor={0,0,127}));
          connect(IC.y, BuckleyLeverett.u1) annotation (points=[-79.6,33; -72.8,
                33; -72.8,26.4; -65.53,26.4],
                                          style(color=74, rgbcolor={0,0,127}));
          connect(BCL.y, BuckleyLeverett.u2) annotation (points=[-79.6,21; -72.8,21;
                -72.8,21.6; -65.53,21.6], style(color=74, rgbcolor={0,0,127}));
        end BuckleyLeverettEquation;
        annotation (Documentation(info="<html>
<p>
This package contains Buckley-Leverett equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end BuckleyLeverett;

      package TransmissionLine
        model TransmissionLinePDE
          Integrator.UniversalIntegrator u(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1) annotation (extent=[-60,20; -20,60]);
          Integrator.UniversalIntegrator v(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1) annotation (extent=[-60,-60; -20,-20]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Transmission line equation</font></h3>
<p>
Implements the Transmission Line equation
</p>

<img align=middle src=\"..\\Images\\tl1.png\">

<p>
where c, h and k are constant values. The initial and boundary conditions are
</p>

<img align=middle src=\"..\\Images\\tl3.png\">

<p>
In order that the integrator block accepts this equation, we must transform
it into two PDEs:
</p>

<img align=middle src=\"..\\Images\\tl2.png\">

<p>
The first equation is implemented in <b>u</b> block, the second in <b>v</b> block.
</p>


</pre>
<p><b>Release Notes: </b></p>
</html>"));
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[32,22; 52,42]);
          Modelica.Blocks.Sources.RealExpression c2[worldModel1.n](y=1.0)
            annotation (extent=[0,16; 20,36]);
          Modelica.Blocks.Math.Add add[worldModel1.n](k2=-1)
            annotation (extent=[64,0; 84,20]);
          Modelica.Blocks.Math.Add add1[worldModel1.n](k2=-1)
            annotation (extent=[60,-38; 80,-18]);
          Modelica.Blocks.Math.Product product1[worldModel1.n]
            annotation (extent=[0,-6; 20,14]);
          Modelica.Blocks.Sources.RealExpression k[worldModel1.n]
            annotation (extent=[-32,-12; -12,8]);
          Modelica.Blocks.Math.Product product2[worldModel1.n]
            annotation (extent=[24,-44; 44,-24]);
          Modelica.Blocks.Sources.RealExpression h[worldModel1.n](y=1.0)
            annotation (extent=[-10,-50; 10,-30]);
          Modelica.Blocks.Sources.RealExpression BCLu
            annotation (extent=[-90,14; -80,28]);
          Modelica.Blocks.Sources.RealExpression BCRu
            annotation (extent=[-90,2; -80,16]);
          Modelica.Blocks.Sources.RealExpression BCLv
            annotation (extent=[-90,-56; -80,-40]);
          Modelica.Blocks.Sources.RealExpression BCRv
            annotation (extent=[-90,-68; -80,-52]);
          Modelica.Blocks.Sources.RealExpression ICv[worldModel1.n]
            annotation (extent=[-90,-44; -80,-28]);
          inner World.worldModel worldModel1
            annotation (extent=[-20,80; 20,100]);
          TLIC tLIC annotation (extent=[-94,30; -78,46]);
          SpaceDerivative.Derivatives.u_xx u_xx
            annotation (extent=[0,42; 20,62]);
        equation
          connect(v.y, u.u) annotation (points=[-18,-28; -6,-28; -6,-12; -98,
                -12; -98,52.4; -61.8,52.4],
                                       style(color=74, rgbcolor={0,0,127}));
          connect(c2.y, product.u2) annotation (points=[21,26; 30,26], style(
                color=74, rgbcolor={0,0,127}));
          connect(u.y, product1.u1) annotation (points=[-18,52; -10,52; -10,10;
                -2,10], style(color=74, rgbcolor={0,0,127}));
          connect(k.y, product1.u2) annotation (points=[-11,-2; -2,-2], style(
                color=74, rgbcolor={0,0,127}));
          connect(product.y, add.u1) annotation (points=[53,32; 56,32; 56,16;
                62,16], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, add.u2) annotation (points=[21,4; 62,4], style(
                color=74, rgbcolor={0,0,127}));
          connect(v.y, product2.u1) annotation (points=[-18,-28; 22,-28], style(
                color=74, rgbcolor={0,0,127}));
          connect(h.y, product2.u2) annotation (points=[11,-40; 22,-40], style(
                color=74, rgbcolor={0,0,127}));
          connect(product2.y, add1.u2) annotation (points=[45,-34; 58,-34],
              style(color=74, rgbcolor={0,0,127}));
          connect(add.y, add1.u1) annotation (points=[85,10; 92,10; 92,-8; 50,
                -8; 50,-22; 58,-22], style(color=74, rgbcolor={0,0,127}));
          connect(add1.y, v.u) annotation (points=[81,-28; 92,-28; 92,-80; -98,
                -80; -98,-27.6; -61.8,-27.6], style(color=74, rgbcolor={0,0,127}));
          connect(BCLu.y, u.u2) annotation (points=[-79.5,21; -72,21; -72,32;
                -61.8,32], style(color=74, rgbcolor={0,0,127}));
          connect(BCRu.y, u.u3) annotation (points=[-79.5,9; -70,9; -70,26;
                -61.8,26], style(color=74, rgbcolor={0,0,127}));
          connect(BCLv.y, v.u2) annotation (points=[-79.5,-48; -61.8,-48],
                                        style(color=74, rgbcolor={0,0,127}));
          connect(BCRv.y, v.u3) annotation (points=[-79.5,-60; -70.75,-60;
                -70.75,-54; -61.8,-54], style(color=74, rgbcolor={0,0,127}));
          connect(ICv.y, v.u1) annotation (points=[-79.5,-36; -70,-36; -70,-42;
                -61.8,-42],             style(color=74, rgbcolor={0,0,127}));
          connect(tLIC.y, u.u1) annotation (points=[-77.2,38; -61.8,38],
                               style(color=74, rgbcolor={0,0,127}));
          connect(u.y, u_xx.u) annotation (points=[-18,52; -2,52], style(color=
                  74, rgbcolor={0,0,127}));
          connect(u_xx.y, product.u1) annotation (points=[21,52; 24,52; 24,38;
                30,38], style(color=74, rgbcolor={0,0,127}));
        end TransmissionLinePDE;

        block TLIC
          extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          inner parameter Integer n = worldModel1.n;

          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          Modelica.Blocks.Sources.Constant const[worldModel1.n](k=1:worldModel1.
                n) annotation (extent=[-80,60; -60,80]);
          Modelica.Blocks.Math.Add add[worldModel1.n](k2=-1)
            annotation (extent=[-40,40; -20,60]);
          annotation (Diagram, Icon(Text(
                extent=[-42,28; 38,-20],
                style(color=3, rgbcolor={0,0,255}),
                string="TLIC")),
            Documentation(info="<html>
<p>
Implements initial condition for the <b>u</b> block of the transmission line equation:
</p>



<img align=middle src=\"..\\Images\\tl4.png\">




</pre>
<p><b>Release Notes: </b></p>
</html>"));
          Modelica.Blocks.Sources.IntegerExpression one[worldModel1.n](y=1)
            annotation (extent=[-80,34; -60,54]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[0,34; 20,54]);
          Modelica.Blocks.Sources.RealExpression pi[worldModel1.n](y=3.14)
            annotation (extent=[-40,16; -20,36]);
          Modelica.Blocks.Math.Division division[worldModel1.n]
            annotation (extent=[34,28; 54,48]);
          Modelica.Blocks.Math.Sin sin[worldModel1.n]
            annotation (extent=[70,-10; 90,10]);
          Modelica.Blocks.Sources.RealExpression realExpression[worldModel1.n](
              y=worldModel1.n - 1) annotation (extent=[0,10; 20,30]);
        equation
          connect(const.y, add.u1) annotation (points=[-59,70; -50,70; -50,56;
                -42,56], style(color=74, rgbcolor={0,0,127}));
          connect(one.y, add.u2) annotation (points=[-59,44; -42,44], style(
                color=45, rgbcolor={255,127,0}));
          connect(add.y, product.u1) annotation (points=[-19,50; -2,50], style(
                color=74, rgbcolor={0,0,127}));
          connect(pi.y, product.u2) annotation (points=[-19,26; -10,26; -10,38;
                -2,38], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, division.u1) annotation (points=[21,44; 32,44],
              style(color=74, rgbcolor={0,0,127}));
          connect(sin.y, y) annotation (points=[91,0; 110,0], style(color=74,
                rgbcolor={0,0,127}));
          connect(division.y, sin.u) annotation (points=[55,38; 60,38; 60,0; 68,
                0], style(color=74, rgbcolor={0,0,127}));
          connect(realExpression.y, division.u2) annotation (points=[21,20; 26,
                20; 26,32; 32,32], style(color=74, rgbcolor={0,0,127}));
        end TLIC;
        annotation (Documentation(info="<html>
<p>
This package contains transmission line equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end TransmissionLine;

      package SimpleSupportedBeam

        model Beam
          Integrator.UniversalIntegrator v(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1) annotation (extent=[-40,20; 0,60]);
          inner World.worldModel worldModel1(u_xx=1, n=60)
            annotation (extent=[-20,88; 20,100]);
          Integrator.UniversalIntegrator w(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1) annotation (extent=[-40,-40; 0,0]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Simple supported beam equation</font></h3>
<p>
Implements the simple supported beam equation
</p>

<img align=middle src=\"..\\Images\\ssb1.png\">

<p>
The initial conditions are
</p>
<p>
<img align=middle src=\"..\\Images\\ssb2.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\ssb3.png\">
</p>

<p>
and boundary conditions are
</p>

<img align=middle src=\"..\\Images\\ssb4.png\">

<p>
Because the integrator block cannot accept the equation in this form, we
transform the PDE above into two first-order PDEs:
</p>

<img align=middle src=\"..\\Images\\ssb5.png\">

<p>
The first equation is implemented in <b>v</b> block, the second in <b>w</b> block. <br>
The analytical solution of this problem is implemented in <b>SSBAnalytic</b> block.
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Gain gain[worldModel1.n](k=-1)
            annotation (extent=[54,-12; 60,-4]);
          Modelica.Blocks.Sources.RealExpression BCLu
            annotation (extent=[-70,26; -60,38]);
          Modelica.Blocks.Sources.RealExpression BCRu
            annotation (extent=[-70,16; -60,28]);
          Modelica.Blocks.Sources.RealExpression ICv[worldModel1.n]
            annotation (extent=[-70,-28; -60,-16]);
          SSBIC sSBIC annotation (extent=[-72,40; -56,48]);
          Modelica.Blocks.Sources.RealExpression BCLv
            annotation (extent=[-70,-38; -60,-26]);
          Modelica.Blocks.Sources.RealExpression BCRv
            annotation (extent=[-70,-48; -60,-36]);
          SSBAnalytic sSBAnalytic annotation (extent=[-40,-80; 0,-60]);
          SpaceDerivative.Derivatives.u_xx u_xx(bcl=0, bcr=0)
            annotation (extent=[20,-18; 40,2]);
          SpaceDerivative.Derivatives.u_xx u_xx1(bcl=-1, bcr=-1)
            annotation (extent=[20,42; 40,62]);
        equation
          connect(BCLu.y,v. u2) annotation (points=[-59.5,32; -41.8,32], style(
                color=74, rgbcolor={0,0,127}));
          connect(BCRu.y,v. u3) annotation (points=[-59.5,22; -52,22; -52,26;
                -41.8,26], style(color=74, rgbcolor={0,0,127}));
          connect(ICv.y,w. u1) annotation (points=[-59.5,-22; -41.8,-22], style(
                color=74, rgbcolor={0,0,127}));
          connect(sSBIC.y,v. u1) annotation (points=[-55.2,44; -48,44; -48,38;
                -41.8,38], style(color=74, rgbcolor={0,0,127}));
          connect(BCLv.y,w. u2) annotation (points=[-59.5,-32; -50,-32; -50,-28;
                -41.8,-28], style(color=74, rgbcolor={0,0,127}));
          connect(BCRv.y,w. u3) annotation (points=[-59.5,-42; -50,-42; -50,-34;
                -41.8,-34], style(color=74, rgbcolor={0,0,127}));
          connect(w.y, u_xx.u) annotation (points=[2,-8; 18,-8], style(color=74,
                rgbcolor={0,0,127}));
          connect(u_xx.y, gain.u) annotation (points=[41,-8; 53.4,-8], style(
                color=74, rgbcolor={0,0,127}));
          connect(gain.y,v. u) annotation (points=[60.3,-8; 68,-8; 68,12; -80,
                12; -80,52.4; -41.8,52.4], style(color=74, rgbcolor={0,0,127}));
          connect(v.y, u_xx1.u) annotation (points=[2,52; 18,52], style(color=
                  74, rgbcolor={0,0,127}));
          connect(u_xx1.y,w. u) annotation (points=[41,52; 68,52; 68,72; -86,72;
                -86,-7.6; -41.8,-7.6], style(color=74, rgbcolor={0,0,127}));
        end Beam;

        block SSBIC
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

        protected
          Real pi = 3.14159265;

        equation
          for i in 1:n loop
            y[i] = sin(pi*(i-1)/(n-1)) + 0.5*(sin(3*pi*(i-1)/(n-1)));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Text(
                extent=[-44,28; 42,-28],
                style(color=3, rgbcolor={0,0,255}),
                string="SSBIC")), Documentation(info="<html>
<p>
Implements the initial condition for the <b>v</b> block of the simple supported beam equation
</p>

<p>
<img align=middle src=\"..\\Images\\ssb2.png\">
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end SSBIC;

        block SSBAnalytic
          extends Icons.BlockIcon4;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

        protected
          Real pi = 3.14159265;

        equation
          for i in 1:n loop
            y[i] = (cos((pi^2)*time))*(sin(pi*(i-1)/(n-1)))  + 0.5*(cos(9*(pi^2)*time))*(sin(3*pi*(i-1)/(n-1)));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Text(
                extent=[-42,26; 64,-42],
                style(color=3, rgbcolor={0,0,255}),
                string="SSBAnalyitc")), Documentation(info="<html>
<p>
Implements the analytical solution of the simple supported beam equation
</p>

<img align=middle src=\"..\\Images\\ssb6.png\">


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end SSBAnalytic;
        annotation (Documentation(info="<html>
<p>
This package contains simple supported beam equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end SimpleSupportedBeam;

      package Poisson
        model PoissonEquation
          inner World.worldModel worldModel1(n=20)
            annotation (extent=[-20,88; 20,100]);
          Integrator.UniversalIntegrator Poisson(
            vb=2,
            ve=worldModel1.n - 1,
            icb=2,
            ice=worldModel1.n - 1,
            bcl=1,
            bcr=1) annotation (extent=[-48,6; 10,64]);
          SpaceDerivative.Derivatives.u_xx u_xx
            annotation (extent=[34,38; 54,58]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Poisson equation</font></h3>
<p>
The poisson problem
</p>

<img align=middle src=\"..\\Images\\poisson1.png\">

<p>
is a time independent problem. In order to be handled from the integrator block, this problem is transformed into a
time dependent problem
</p>

<img align=middle src=\"..\\Images\\poisson2.png\">

<p>
with initial condition
</p>

<p>
<img align=middle src=\"..\\Images\\poisson3.png\">
</p>

<p>
The solution of the time dependent problem approaches with time the solution of the time independent problem. <br>
The analytical solution of this problem is implemented in <b>PAN</b> block.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Add add[worldModel1.n](k1=+1, k2=+1)
            annotation (extent=[68,22; 82,36]);
          Modelica.Blocks.Sources.RealExpression BCL
            annotation (extent=[-86,16; -70,32]);
          Modelica.Blocks.Sources.RealExpression BCR
            annotation (extent=[-86,4; -70,20]);
          PIC pIC annotation (extent=[-88,34; -68,54]);
          PAN pAN annotation (extent=[-40,-80; 0,-60]);
          Modelica.Blocks.Sources.Constant const[worldModel1.n](k=1:worldModel1.
                 n) annotation (extent=[-48,-12; -36,0]);
          Modelica.Blocks.Sources.IntegerExpression one[worldModel1.n](y=1)
            annotation (extent=[-48,-34; -36,-18]);
          Modelica.Blocks.Math.Add add1[worldModel1.n](k2=-1)
            annotation (extent=[-24,-22; -12,-10]);
          Modelica.Blocks.Math.Product product[worldModel1.n]
            annotation (extent=[-2,-30; 10,-18]);
          Modelica.Blocks.Sources.RealExpression pi[worldModel1.n](y=3.14)
            annotation (extent=[-24,-40; -12,-24]);
          Modelica.Blocks.Math.Division division[worldModel1.n]
            annotation (extent=[20,-38; 32,-26]);
          Modelica.Blocks.Sources.RealExpression deltax[worldModel1.n](y=
                worldModel1.n - 1) annotation (extent=[-2,-48; 10,-32]);
          Modelica.Blocks.Math.Sin sin[worldModel1.n]
            annotation (extent=[40,-38; 52,-26]);
        equation
          connect(Poisson.y, u_xx.u) annotation (points=[12.9,52.4; 22.45,52.4;
                22.45,48; 32,48], style(color=74, rgbcolor={0,0,127}));
          connect(u_xx.y, add.u1) annotation (points=[55,48; 60,48; 60,33.2;
                66.6,33.2], style(color=74, rgbcolor={0,0,127}));
          connect(add.y, Poisson.u) annotation (points=[82.7,29; 90,29; 90,70;
                -62,70; -62,52.98; -50.61,52.98], style(color=74, rgbcolor={0,0,
                  127}));
          connect(BCR.y, Poisson.u3) annotation (points=[-69.2,12; -60,12; -60,
                14.7; -50.61,14.7], style(color=74, rgbcolor={0,0,127}));
          connect(pIC.y, Poisson.u1) annotation (points=[-67,44; -60,44; -60,
                32.1; -50.61,32.1], style(color=74, rgbcolor={0,0,127}));
          connect(const.y, add1.u1) annotation (points=[-35.4,-6; -30,-6; -30,
                -12.4; -25.2,-12.4], style(color=74, rgbcolor={0,0,127}));
          connect(one.y, add1.u2) annotation (points=[-35.4,-26; -30,-26; -30,
                -19.6; -25.2,-19.6], style(color=45, rgbcolor={255,127,0}));
          connect(add1.y, product.u1) annotation (points=[-11.4,-16; -6,-16; -6,
                -20.4; -3.2,-20.4], style(color=74, rgbcolor={0,0,127}));
          connect(pi.y, product.u2) annotation (points=[-11.4,-32; -6,-32; -6,
                -27.6; -3.2,-27.6], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, division.u1) annotation (points=[10.6,-24; 14,-24;
                14,-28.4; 18.8,-28.4], style(color=74, rgbcolor={0,0,127}));
          connect(deltax.y, division.u2) annotation (points=[10.6,-40; 14,-40;
                14,-35.6; 18.8,-35.6], style(color=74, rgbcolor={0,0,127}));
          connect(division.y, sin.u) annotation (points=[32.6,-32; 38.8,-32],
              style(color=74, rgbcolor={0,0,127}));
          connect(sin.y, add.u2) annotation (points=[52.6,-32; 60,-32; 60,24.8;
                66.6,24.8], style(color=74, rgbcolor={0,0,127}));
          connect(BCL.y, Poisson.u2) annotation (points=[-69.2,24; -59.905,24;
                -59.905,23.4; -50.61,23.4], style(color=74, rgbcolor={0,0,127}));
        end PoissonEquation;

        block PIC
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

        equation
          for i in 1:n loop
            y[i] = sin(3*3.14*((i-1)/n));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Text(
                extent=[-52,28; 50,-26],
                style(color=3, rgbcolor={0,0,255}),
                string="PIC")), Documentation(info="<html>
<p>
Implements the initial condition for the time dependent poisson equation
</p>

<p>
<img align=middle src=\"..\\Images\\poisson3.png\">
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end PIC;

        block PAN
          extends Icons.BlockIcon4;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;

        equation
          for i in 1:n loop
            y[i] = (1/(3.14)^2)*sin(3.14*((i-1)/n));
          end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Text(
                extent=[-56,34; 58,-34],
                style(color=3, rgbcolor={0,0,255}),
                string="PAN")), Documentation(info="<html>
<p>
Implements the analytical solution of the poisson equation
</p>

<img align=middle src=\"..\\Images\\poisson4.png\">


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end PAN;
        annotation (Documentation(info="<html>
<p>
This package contains Poisson equation solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Poisson;
      annotation (Documentation(info="<html>
<p>
This package contains examples of partial differential equations solved with the Method of Lines.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end Examples;
    annotation (Documentation(info="<html>
<p>
This package contains necessary blocks for solving partial differential equations with Method of Lines.
To understand the use of the blocks, many examples are implemented (PDE->MOL->Examples).
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
  end MOL;

  package FiniteVolume
    package FVMIntegrator

      block FVIntegrator

        extends Icons.BlockIcon1;

      outer PDE.World.worldModel worldModel1;
      parameter Integer n = worldModel1.n;
      //parameter Real beta = 1;
      parameter Integer vb = gcl+1 "|Unknowns| The left most unknown";
      parameter Integer ve = gcl + worldModel1.n
          "|Unknowns| The right most unknown";

      parameter Integer icb = gcl+1
          "|Initial Condition| Begin of the initial condition";
      parameter Integer ice = gcl + worldModel1.n
          "|Initial Condition| End of the initial condition";

      parameter Integer bcl = 1
          "|Boundary Conditions| Boundary condition at the left (0: no; 1: yes)";
      parameter Integer bcr = 1
          "|Boundary Conditions| Boundary condition at the right (0: no; 1: yes)";

      parameter Integer gcl = 2
          "|Boundary Conditions| Number of ghost cells at the left";
      parameter Integer gcr = 2
          "|Boundary Conditions| Number of ghost cells at the right";

      parameter Real delta_x = 1/n;
      //parameter Real delta_t = 0.1;

      Real q[n+gcl+gcr];

      equation
        y = q;

        for i in 1:gcl loop
           q[i]= u3[i];
        end for;

        if bcl == 1 then
           q[gcl+1] = q[gcl];
        end if;

        for i in 1:gcr loop
           q[gcl+n+i] = u4[i];
        end for;

        if bcr == 1 then
           q[gcl+n] = q[gcl+n+1];
        end if;

        for i in vb:ve loop
          der(q[i]) = -(1/delta_x)*(u[i-gcl+1]-u[i-gcl]);
        end for;

      initial equation

        for i in icb:ice loop
          q[i] =  u2[i-gcl];
        end for;

      // equation
      //   y = q;
      //
      //   if bcl == 1 then
      //     for i in 1:gcl loop
      //       q[i]= u3[i];
      //     end for;
      //   end if;
      //
      //   if bcr == 1 then
      //     for i in 1:gcr loop
      //       q[gcl+n+i] = u4[i];
      //     end for;
      //   end if;
      //
      //   for i in vb:ve loop
      //     der(q[i]) = -(1/delta_x)*(u[i-gcl+1]-u[i-gcl]);
      //   end for;
      //
      // initial equation
      //
      //   for i in icb:ice loop
      //     q[i] =  u2[i-gcl];
      //   end for;

      public
        Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
          annotation (extent=[-122,66; -100,94]);
        Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + worldModel1.gcl
           + worldModel1.gcr]
          annotation (extent=[100,50; 120,70]);
        Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n]
          annotation (extent=[-122,-14; -100,14]);
        Modelica.Blocks.Interfaces.RealInput u3[worldModel1.gcl]
          annotation (extent=[-122,-64; -100,-36]);
        Modelica.Blocks.Interfaces.RealInput u4[worldModel1.gcr]
          annotation (extent=[-122,-94; -100,-66]);
        annotation (Diagram, Icon(
            Text(
              extent=[-40,162; 40,124],
              style(color=3, rgbcolor={0,0,255}),
              string="%name"),
            Text(
              extent=[-108,88; -60,70],
              style(color=3, rgbcolor={0,0,255}),
              string="F"),
            Text(
              extent=[-100,10; -64,-10],
              style(color=3, rgbcolor={0,0,255}),
              string="IC"),
            Text(
              extent=[-102,-38; -58,-58],
              style(color=3, rgbcolor={0,0,255}),
              string="gcl"),
            Text(
              extent=[-108,-68; -50,-88],
              style(color=3, rgbcolor={0,0,255}),
              string="gcr"),
            Text(
              extent=[56,74; 90,46],
              style(color=3, rgbcolor={0,0,255}),
              string="Q")),
          Documentation(info="<html>
<p>
Implements the cell average update rule. In one dimension, the finite volume method subdivide the domain into <br>
cells (intervals) and approximates the integral of the unknown function <b>q</b> over each of these cells at each time step (see figure below). <br>
The <b>ghost cells</b> are the boundary cells that are introduced to avoid writing special formulas for the boundary cells.
</p>

<img align=middle src=\"..\\Images\\fvmi7.png\">

<p>
Denote the i-th cell by
</p>

<img align=middle src=\"..\\Images\\fvmi.png\">

<p>
Then the approximation to the average of <b>q</b> in the cell C<sub>i</sub> at time t, which we denote with Q<sub>i</sub> is
</p>

<img align=middle src=\"..\\Images\\fvmi1.png\">


<p>
The approximation to this average derives from the integral form of the conservation law
</p>

<img align=middle src=\"..\\Images\\fvmi2.png\">

<p>
which states that the average within the cell can only changes due to the fluxes at the boundaries (if we assume that <br>
no source or sink is present in the cell). <br>
If we integrate this expression in time from <i>t</i> to <i>t+deltat</i>, we obtain <br>
</p>

<img align=middle src=\"..\\Images\\fvmi3.png\">

<p>
and dividing by <i>deltax</i> we reach the form
</p>

<img align=middle src=\"..\\Images\\fvmi4.png\">

<p>
which give us an explicit time marching algorithm. By using the notation for averages introduced above we can write
</p>

<img align=middle src=\"..\\Images\\fvmi5.png\">

<p>
where
</p>

<img align=middle src=\"..\\Images\\fvmi6.png\">

<p>
approximates the average flux along the interface x<sub>i-1/2</sub>. <br>
The FVMIntegrator block implements this average update rule. Initial condition of the problem can be passed to the <i>IC</i> input of the block, <br>
whereas the boundary conditions to the corresponding ghost cells. The number of ghost cells depends on the method we use. In the present <br>
package, two ghost cells at the left and at the right are enough for every method implemented.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end FVIntegrator;

      annotation (Documentation(info="<html>
<p>
This package contains integrator block that implements Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end FVMIntegrator;

    package Fluxes

      package DiffusionFlux
        block DiffusionFlux
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Real beta = 1.0;
        parameter Real deltax = 1/n;

        equation
          // for i in 1:gcl loop
          //   y[i] = u[i];
          // end for;

          for i in gcl:gcl+n loop
            y[i-gcl+1] = (-beta/(deltax))*(u[i+1] - u[i]);
          end for;

          // for i in gcl+n+1:gcl+n+gcr loop
          //   y[i] = u[i];
          // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.
            gcl + worldModel1.gcr] annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-38,46; 42,-40],
                style(color=3, rgbcolor={0,0,255}),
                string="Diffusion Flux")),
            Documentation(info="<html>
<p>
Implements the diffusion flux
</p>

<img align=middle src=\"..\\Images\\df1.png\">

<p>
By using this flux, the average update rule becomes:
</p>

<img align=middle src=\"..\\Images\\df2.png\">



</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end DiffusionFlux;
        annotation (Documentation(info="<html>
<p>
This package contains the DiffusionFlux block for the computation of fluxes at the cell boundaries.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end DiffusionFlux;

      package UpwindFlux
        block Upwind
          extends Icons.BlockIcon5;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;

        equation
        if u1 > 0 then

          for i in 1:n+1 loop
            y[i] = u1*u[gcl+i-1];
          end for;

        else

          for i in 1:n+1 loop
            y[i] = u1*u[gcl+i];
          end for;

        end if;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcr
             + worldModel1.gcr] annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Documentation(info="<html>
<p>
Implements the upwind flux: F<sub>i-1/2</sub> = uQ<sub>i-1</sub> if u > 0 and F<sub>i-1/2</sub> = uQ<sub>i</sub> otherwise. <br>
By using this flux the average update rule becomes
</p>

<img align=middle src=\"..\\Images\\uw1.png\">

<p>
if u > 0, and
</p>

<img align=middle src=\"..\\Images\\uw2.png\">

<p>
otherwise.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"),  Icon(Text(
                extent=[-44,38; 44,-34],
                style(color=3, rgbcolor={0,0,255}),
                string="UW")));
          Modelica.Blocks.Interfaces.RealInput u1
            annotation (extent=[-140,-80; -100,-40]);
        end Upwind;
        annotation (Documentation(info="<html>
<p>
This package contains the Upwind block for the computation of fluxes at the cell boundaries.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end UpwindFlux;

      package LaxFriedrichFlux
        block LF
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        // parameter Integer gcl = worldModel1.gcl;
        // parameter Integer gcr = worldModel1.gcr;
        parameter Real alpha = 0.0;

        equation
        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;

        for i in 1:n+1 loop
          y[i] = 0.5*(u[i] + u1[i]) - 0.5*alpha*(u2[i] - u3[i]);
        end for;

        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
                                annotation (extent=[-118,68; -100,92]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
                                annotation (extent=[-118,38; -100,62]);
          Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
                                annotation (extent=[-118,-62; -100,-38]);
          Modelica.Blocks.Interfaces.RealInput u3[worldModel1.n + 1]
                                annotation (extent=[-118,-92; -100,-68]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
                                annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-36,36; 36,-30],
                style(color=3, rgbcolor={0,0,255}),
                string="LF")),
            Documentation(info="<html>

<p>
Implements the Lax-Friedrichs flux
</p>

<img align=middle src=\"..\\Images\\lf1.png\">

<p>
By using this flux, the average update rule becomes:
</p>

<img align=middle src=\"..\\Images\\lf2.png\">



</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end LF;
        annotation (Documentation(info="<html>
<p>
This package contains the LF block for the computation of fluxes at the cell boundaries.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end LaxFriedrichFlux;

      package LaxWendroffFlux

        block LaxWendroff
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Real deltax = 1/n;
        parameter Real deltat = worldModel1.deltat;

        equation
          for i in 1:n+1 loop
            y[i] = 0.5*u1*(u[gcl+i] + u[gcl+i-1]) - 0.5*(u1^2)*(deltat/deltax)*(u[gcl+i] - u[gcl+i-1]);
          end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcr
             + worldModel1.gcr] annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-42,46; 44,-38],
                style(color=3, rgbcolor={0,0,255}),
                string="LW")),
            Documentation(info="<html>
<p>
Implements the Lax-Wendroff flux. By using this flux, the average update rule becomes:
</p>

<img align=middle src=\"..\\Images\\lw.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end LaxWendroff;
        annotation (Documentation(info="<html>
<p>
This package contains the LaxWendroffFlux block for the computation of fluxes at the cell boundaries.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end LaxWendroffFlux;

      package Roe
        package Averages
          block Daverage
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;

          equation
          for i in 1:n+1 loop
            y[i] = sqrt(u[i]*u1[i]);
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-140,-80; -100,-40]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            annotation (Diagram, Icon(Bitmap(extent=[-58,58; 54,-46], name=
                      "Images/rtilde.png")),
              Documentation(info="<html>
<p>
Takes the reconstructed values rho<sub>i</sub><sup>+</sup> and rho<sub>i</sub><sup>-</sup> and computes the average
</p>

<img align=middle src=\"..\\Images\\r2.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          end Daverage;

          block Vaverage
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;

            Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
              annotation (extent=[-140,-60; -100,-20]);
            Modelica.Blocks.Interfaces.RealInput u3[worldModel1.n + 1]
              annotation (extent=[-140,-100; -100,-60]);
          equation
          for i in 1:n+1 loop
            y[i] = ((sqrt(u[i]))*u1[i] + (sqrt(u2[i]))*u3[i])/((sqrt(u[i])) + (sqrt(u2[i])));
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-140,60; -100,100]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-140,20; -100,60]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            annotation (Diagram, Icon(Bitmap(extent=[-62,70; 62,-40], name=
                      "Images/vtilde.png")),
              Documentation(info="<html>
<p>
Takes the reconstructed values rho<sub>i</sub><sup>+</sup>, rho<sub>i</sub><sup>-</sup>, v<sub>i</sub><sup>+</sup> and v<sub>i</sub><sup>-</sup> and computes the average
</p>

<img align=middle src=\"..\\Images\\r3.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          end Vaverage;

          block Haverage
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;

          equation
          for i in 1:n+1 loop
            y[i] = ((sqrt(u[i]))*u1[i] + (sqrt(u2[i]))*u3[i])/((sqrt(u[i])) + (sqrt(u2[i])));
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-140,60; -100,100]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-140,20; -100,60]);
            Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
              annotation (extent=[-140,-60; -100,-20]);
            Modelica.Blocks.Interfaces.RealInput u3[worldModel1.n + 1]
              annotation (extent=[-140,-100; -100,-60]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            annotation (Diagram, Icon(Bitmap(extent=[-68,74; 64,-42], name=
                      "Images/htilde.png")),
              Documentation(info="<html>
<p>
Takes the reconstructed values rho<sub>i</sub><sup>+</sup>, rho<sub>i</sub><sup>-</sup>, h<sub>i</sub><sup>+</sup> and hv<sub>i</sub><sup>-</sup> and computes the average
</p>

<img align=middle src=\"..\\Images\\r4.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          end Haverage;

          block Aaverage
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Real lambda = 1.4;

          equation
          for i in 1:n+1 loop
            y[i] = sqrt((lambda - 1)*(u[i] - 0.5*((u1[i])^2)));
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-140,-80; -100,-40]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            annotation (Diagram, Icon(Bitmap(extent=[-72,68; 64,-56], name=
                      "Images/atilde.png")),
              Documentation(info="<html>
<p>
Takes the averages h and v and computes the average
</p>

<img align=middle src=\"..\\Images\\r5.png\">

<p>
where gamma = 1.4.
</p>
</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          end Aaverage;
          annotation (Documentation(info="<html>
<p>
This package contains blocks for the computation of roe큦 averages
</p>

<img align=middle src=\"..\\Images\\r1.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end Averages;

        package DeltaU
          block Deltau
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;

          equation
          for i in 1:n+1 loop
            y[i] = -(u[i] - u1[i]);
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-124,58; -100,86]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-124,-86; -100,-58]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            annotation (Diagram, Icon(Bitmap(extent=[-76,74; 74,-64], name=
                      "Images/deltau.png")),
              Documentation(info="<html>
<p>
Takes the reconstructed values + and - and computes their difference. For example for pressure p we have
</p>
<img align=middle src=\"..\\Images\\r7.png\">
</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          end Deltau;
          annotation (Documentation(info="<html>
<p>
This package contains DeltaU block that computes the difference between the reconstructed values + and -.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end DeltaU;

        package Lambda
          block Lambdas
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          inner parameter Integer n = worldModel1.n;

            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-124,44; -100,78]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-124,-78; -100,-44]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,50; 120,70]);
            Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            Modelica.Blocks.Interfaces.RealOutput y2[worldModel1.n + 1]
              annotation (extent=[100,-70; 120,-50]);
            annotation (Diagram, Icon(
                Bitmap(extent=[-76,50; 68,-44], name="Images/lambda.png"),
                Bitmap(extent=[-100,76; -62,50], name="Images/vtilde.png"),
                Bitmap(extent=[-96,-42; -62,-74], name="Images/atilde.png"),
                Bitmap(extent=[56,76; 90,44], name="Images/lambda1.png"),
                Bitmap(extent=[52,14; 92,-14], name="Images/lambda2.png"),
                Bitmap(extent=[52,-48; 90,-74], name="Images/lambda3.png")),
              Documentation(info="<html>
<p>
Takes the averages v, a and h and computes the three eigenvalues
</p>
<img align=middle src=\"..\\Images\\r6.png\">
</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
            Modelica.Blocks.Math.Add add[worldModel1.n + 1](k2=-1)
              annotation (extent=[-56,44; -36,64]);
            Modelica.Blocks.Math.Add add1[worldModel1.n + 1]
              annotation (extent=[-60,-64; -40,-44]);
          equation
            connect(u, add.u1) annotation (points=[-112,61; -85,61; -85,60; -58,
                  60],                                               style(color=
                    74, rgbcolor={0,0,127}));
            connect(u, add1.u1) annotation (points=[-112,61; -80,61; -80,-48;
                  -62,-48],
                        style(color=74, rgbcolor={0,0,127}));
            connect(u, y1) annotation (points=[-112,61; -80,61; -80,0; 110,0],
                style(color=74, rgbcolor={0,0,127}));
            connect(u1, add.u2) annotation (points=[-112,-61; -88,-61; -88,48;
                  -58,48], style(color=74, rgbcolor={0,0,127}));
            connect(u1, add1.u2) annotation (points=[-112,-61; -87,-61; -87,-60;
                  -62,-60],                                              style(
                  color=74, rgbcolor={0,0,127}));
            connect(add.y, y) annotation (points=[-35,54; 80,54; 80,60; 110,60],
                style(color=74, rgbcolor={0,0,127}));
            connect(add1.y, y2) annotation (points=[-39,-54; 80,-54; 80,-60; 110,
                  -60], style(color=74, rgbcolor={0,0,127}));
          end Lambdas;
          annotation (Documentation(info="<html>
<p>
This package contains Lambdas block that computes the eigenvalues of the Euler system.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end Lambda;

        package Wave
          block Waves
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;

          equation
          for j in 1:n+1 loop
            y[1, j] = 1;
            y[2, j] = u[j] - u1[j];
            y[3, j] = u2[j] - u[j]*u1[j];

            y1y[1, j] = 1;
            y1y[2, j] = u[j];
            y1y[3, j] = 0.5*(u[j])^2;

            y2y[1, j] = 1;
            y2y[2, j] = u[j] + u1[j];
            y2y[3, j] = u2[j] + u[j]*u1[j];
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-126,48; -100,74]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-126,-12; -100,14]);
            Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
              annotation (extent=[-126,-74; -100,-46]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,50; 120,70]);
            Modelica.Blocks.Interfaces.RealOutput y1y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-10; 120,10]);
            Modelica.Blocks.Interfaces.RealOutput y2y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-70; 120,-50]);
            annotation (Diagram, Icon(
                Text(
                  extent=[58,70; 92,50],
                  string="W1",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[62,10; 86,-10],
                  string="W2",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[62,-44; 86,-78],
                  string="W3",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[-30,28; 30,-28],
                  string="W",
                  style(color=0, rgbcolor={0,0,0})),
                Bitmap(extent=[-100,78; -50,48], name="Images/vtilde.png"),
                Bitmap(extent=[-94,20; -56,-18], name="Images/atilde.png"),
                Bitmap(extent=[-100,-42; -48,-74], name="Images/htilde.png")),
              Documentation(info="<html>
<p>
Takes the averages v, a and h and computes the eigenvectors
</p>
<img align=middle src=\"..\\Images\\r8.png\">
<img align=middle src=\"..\\Images\\r9.png\">
<img align=middle src=\"..\\Images\\r10.png\">
</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          end Waves;
          annotation (Documentation(info="<html>
<p>
This package contains Waves block that computes the eigenvectors of the Euler system.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end Wave;

        package WaveStrength
          block a

            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;

          equation
          for j in 1:n+1 loop
            //Inserting 0.0001 to avoid division by zero
            y[j] = (1/(2*((u1[j]+0.0001)^2)))*(u4[j] - u[j]*u1[j]*u3[j]);
            y1[j] = (1/(((u1[j]+0.0001)^2)))*((u1[j]^2)*u2[j] - u4[j]);
            y2[j] = (1/(2*((u1[j]+0.0001)^2)))*(u4[j] + u[j]*u1[j]*u3[j]);
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-118,66; -100,94]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-118,24; -100,54]);
            Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
              annotation (extent=[-118,-14; -100,16]);
            Modelica.Blocks.Interfaces.RealInput u3[worldModel1.n + 1]
              annotation (extent=[-118,-54; -100,-26]);
            Modelica.Blocks.Interfaces.RealInput u4[worldModel1.n + 1]
              annotation (extent=[-118,-92; -100,-64]);
            annotation (Icon(
                Text(
                  extent=[70,68; 98,54],
                  string="a1",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[66,8; 102,-6],
                  string="a2",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[72,-54; 96,-68],
                  string="a3",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[-34,32; 42,-24],
                  string="a",
                  style(color=0, rgbcolor={0,0,0})),
                Bitmap(extent=[-100,98; -54,64], name="Images/rtilde.png"),
                Bitmap(extent=[-98,62; -56,20], name="Images/atilde.png"),
                Bitmap(extent=[-96,22; -48,-24], name="Images/deltar.png"),
                Bitmap(extent=[-100,-22; -46,-60], name="Images/deltav.png"),
                Bitmap(extent=[-100,-62; -46,-98], name="Images/deltap.png")),
                                     Diagram,
              Documentation(info="<html>
<p>
Takes the averages rho, a and differences deltarho, deltav and deltap as input and computes the three wave strengths
</p>

<img align=middle src=\"..\\Images\\r11.png\">

</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,50; 120,70]);
            Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            Modelica.Blocks.Interfaces.RealOutput y2[worldModel1.n + 1]
              annotation (extent=[100,-70; 120,-50]);
          end a;
          annotation (Documentation(info="<html>
<p>
This package contains <i>a</i> block that computes the wave strengths.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end WaveStrength;

        package FluxDifference
          block FluxDiff
          extends Icons.BlockIcon;

            outer PDE.World.worldModel worldModel1;
            parameter Integer n=worldModel1.n;
            parameter Integer m=worldModel1.m;

          equation
          for j in 1:n+1 loop
          if u2[j] > 0 then
            for i in 1:m loop
              y[i, j] = u[j]*u1[i, j];
              y1[i, j] = 0.0;
            end for;
          else
            for i in 1:m loop
              y1[i, j] = u[j]*u1[i, j];
              y[i, j] = 0.0;
            end for;
          end if;
          end for;

            // for j in 1:n+1 loop
            //   for i in 1:m loop
            //     if u2[j] > 0 then
            //       y[i, j] = u[j]*u1[i, j];
            //       y1[i, j] = 0.0;
            //     else
            //       y1[i, j] = u[j]*u1[i, j];
            //       y[i, j] = 0.0;
            //     end if;
            //   end for;
            // end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-118,50; -100,76]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[-118,-12; -100,14]);
            Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
              annotation (extent=[-118,-72; -100,-44]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,50; 120,70]);
            Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-70; 120,-50]);
            annotation (Diagram, Icon(
                Text(
                  extent=[-90,74; -54,50],
                  string="a",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[-92,14; -48,-14],
                  string="W",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[50,74; 90,46],
                  string="W+",
                  style(color=0, rgbcolor={0,0,0})),
                Text(
                  extent=[54,-38; 84,-80],
                  string="W-",
                  style(color=0, rgbcolor={0,0,0})),
                Bitmap(extent=[-102,-40; -38,-80], name="Images/lambda.png")),
              Documentation(info="<html>
<p>
Computes the waves
</p>

<img align=middle src=\"..\\Images\\r12.png\">

<p>
and
</p>

<img align=middle src=\"..\\Images\\r13.png\">

<p>
according to the eigenvalue lambda.
</p>
</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          end FluxDiff;
          annotation (Documentation(info="<html>
<p>
This package contains FluxDiff block that computes the waves W+ and W-.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end FluxDifference;

        package IntegratorRoe
          block Integrator
            extends Icons.BlockIcon1;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer m = worldModel1.m;

          parameter Integer vb = gcl+1 "|Unknowns| The left most unknown";
          parameter Integer ve = gcl + worldModel1.n
              "|Unknowns| The right most unknown";

          parameter Integer icb = gcl+1
              "|Initial Condition| Begin of the initial condition";
          parameter Integer ice = gcl + worldModel1.n
              "|Initial Condition| End of the initial condition";

          parameter Integer bcl = 1
              "|Boundary Conditions| Boundary condition at the left (0: no; 1: yes)";
          parameter Integer bcr = 1
              "|Boundary Conditions| Boundary condition at the right (0: no; 1: yes)";

          parameter Integer gcl = 2
              "|Boundary Conditions| Number of ghost cells at the left";
          parameter Integer gcr = 2
              "|Boundary Conditions| Number of ghost cells at the right";

          parameter Real delta_x = 1/n;
          //parameter Real deltat = 0.00001;
          parameter Integer index = 1;

          Real q[m, n + gcl + gcr];

            Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
              1] annotation (extent=[-110,68; -100,88]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n]
              annotation (extent=[-110,-18; -100,2]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + worldModel1.gcl + worldModel1.gcr]
                                     annotation (extent=[100,50; 120,70]);
            Modelica.Blocks.Interfaces.RealInput u2[worldModel1.m,worldModel1.gcl]
              annotation (extent=[-110,-62; -100,-42]);
            Modelica.Blocks.Interfaces.RealInput u3[worldModel1.m,worldModel1.gcr]
              annotation (extent=[-110,-88; -100,-68]);
            Modelica.Blocks.Interfaces.RealInput u4[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[-110,40; -100,60]);
          equation
            y = q;

          for i in 1:m loop
            for j in 1:gcl loop
               q[i, j]= u2[i, j];
            end for;
          end for;

          for i in 1:m loop
            if bcl == 1 then
               q[i, gcl+1] = q[i, gcl];
            end if;
          end for;

          for i in 1:m loop
            for j in 1:gcr loop
               q[i, gcl+n+j] = u3[i, j];
            end for;
          end for;

          for i in 1:m loop
            if bcr == 1 then
               q[i, gcl+n] = q[i, gcl+n+1];
            end if;
          end for;

          for j in vb:ve loop
            for i in 1:m loop
              der(q[i, j]) = -(1/delta_x)*(u[i, j-gcl] + u4[i, j-gcl+1]);
            end for;
          end for;

          // for j in 1:m loop
          //   for i in vb:ve loop
          //     der(q[j, i]) = -(1/delta_x)*(u[j, i-gcl] + u4[j, i-gcl+1]);
          //   end for;
          // end for;

          // when sample(0, deltat) then
          // for j in 1:m loop
          //   for i in vb:ve loop
          //     q[j, i] = pre(q[j, i]) -(deltat/delta_x)*(pre(u[j, i-gcl]) + pre(u4[j, i-gcl+1]));
          //     //der(q[j, i]) = -(1/delta_x)*(u[j, i-gcl] + u4[j, i-gcl+1]);
          //   end for;
          // end for;
          // end when;

          initial equation

          for j in 1:m loop
            for i in icb:ice loop
              q[j, i] =  u1[j, i-gcl];
            end for;
          end for;

            annotation (Diagram, Icon(Text(
                  extent=[-48,48; 48,-42],
                  style(color=3, rgbcolor={0,0,255}),
                  string="%name")),
              Documentation(info="<html>
<p>
The Flux Limiter Integrator takes as input the fluctuation matrices <b>A<sup>+</sup></b> and <b>A<sup>-</sup></b>, initial condition matrix <b>IC</b> and two boundary condition matrices <b>gcl</b> and <b>gcr</b>. <br>
As output it computes the cell averages at the next time step by using the formula
</p>

<img align=middle src=\"..\\Images\\r14.png\">

<p>
where
</p>

<img align=middle src=\"..\\Images\\r15.png\">

<p>
are the sum of all right-going waves, respectively left-going waves, from the interfaces x<sub>i-1/2</sub>.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          end Integrator;
          annotation (Documentation(info="<html>
<p>
This package contains integrator block that implements Finite Volume Methods with Roe큦 flux.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end IntegratorRoe;

        package AverageFilter
          block Filter
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;
          parameter Integer row = 1;

          equation
          for j in 1:n+gcl+gcr loop
              y[j] = u[row, j];
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
              worldModel1.gcl + worldModel1.gcr]
              annotation (extent=[-140,-20; -100,20]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + worldModel1.
              gcl + worldModel1.gcr] annotation (extent=[100,-10; 120,10]);
            annotation (Diagram, Icon(Text(
                  extent=[-36,26; 38,-18],
                  style(color=3, rgbcolor={0,0,255}),
                  string="%name")),
              Documentation(info="<html>
<p>
This package contains Filter block that filters the requested average vector from the average matrix.
This is done by specifying the value of the <i>row</i> parameter.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          end Filter;
          annotation (Documentation(info="<html>
<p>
This package contains Filter block that filters the requested average vector from the average matrix.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end AverageFilter;

        annotation (Documentation(info="<html>
<p>
This package contains blocks for the computation of Roe큦 flux.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Roe;
      annotation (Documentation(info="<html>
<p>
This package contains different flux blocks.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end Fluxes;

    package LDLR
      package d
        block D1minus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Real deltax = 1/n;

        equation
        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;

        for i in gcl:gcl+n loop
          y[i-gcl+1] = (u[i] - u[i-1])/deltax;
        end for;

        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-28,20; 38,-16],
                style(color=3, rgbcolor={0,0,255}),
                string="D1- ")),
            Documentation(info="<html>
<p>
Implements the lateral derivative d<sub>1</sub> for u<sup>-</sup>(x<sub>i+1/2</sub>)
</p>

<img align=middle src=\"..\\Images\\log5.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end D1minus;

        block D2minus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Real deltax = 1/n;

        equation
        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;

        for i in gcl:gcl+n loop
          y[i-gcl+1] = (u[i+1] - u[i])/deltax;
        end for;

        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-26,20; 38,-16],
                style(color=3, rgbcolor={0,0,255}),
                string="D2- ")),
            Documentation(info="<html>
<p>
Implements the lateral derivative d<sub>2</sub> for u<sup>-</sup>(x<sub>i+1/2</sub>)
</p>

<img align=middle src=\"..\\Images\\log6.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end D2minus;

        block D1plus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Real deltax = 1/n;

        equation
        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;

        for i in gcl:gcl+n loop
          y[i-gcl+1] = (u[i+1] - u[i])/deltax;
        end for;

        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-30,18; 42,-14],
                style(color=3, rgbcolor={0,0,255}),
                string="D1+ ")),
            Documentation(info="<html>
<p>
Implements the lateral derivative d<sub>1</sub> for u<sup>+</sup>(x<sub>i+1/2</sub>)
</p>

<img align=middle src=\"..\\Images\\log7.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end D1plus;

        block D2plus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Real deltax = 1/n;

        equation
        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;

        for i in gcl:gcl+n loop
          y[i-gcl+1] = (u[i+2] - u[i+1])/deltax;
        end for;

        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-26,18; 40,-16],
                style(color=3, rgbcolor={0,0,255}),
                string="D2+ ")),
            Documentation(info="<html>
<p>
Implements the lateral derivative d<sub>2</sub> for u<sup>+</sup>(x<sub>i+1/2</sub>)
</p>

<img align=middle src=\"..\\Images\\log8.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end D2plus;

        annotation (Documentation(info="<html>
<p>
Implements the lateral derivatives d<sub>1</sub> and d<sub>2</sub>
</p>

<img align=middle src=\"..\\Images\\log4.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end d;

      package c
        block C1

          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        // parameter Integer gcl = worldModel1.gcl;
        // parameter Integer gcr = worldModel1.gcr;
        parameter Real deltax = 1/n;

        parameter Real q = 1.4;
        parameter Real tol = 0.1*(deltax^q);

        protected
        Real u_[worldModel1.n + 1];
        Real u1_[worldModel1.n + 1];

        equation
        for i in 1:n+1 loop
          u_[i] = noEvent(if u[i] < 0 then -u[i] else u[i]);
          u1_[i] = noEvent(if u1[i] < 0 then -u1[i] else u1[i]);

         y[i] = (1-tol)*(1 + tol - ((2*((u_[i])^q)*((u1_[i])^q) + tol)/((u_[i])^(2*q) + (u1_[i])^(2*q) + tol)));

         //y[i] = (1-tol)*(1 + tol - ((2*((sqrt(u[i]*u[i]))^q)*((sqrt(u1[i]*u1[i]))^q) + tol)/((sqrt(u[i]*u[i]))^(2*q) + (sqrt(u1[i]*u1[i]))^(2*q) + tol)));

        // y[i] = (1-tol)*(1 + tol - ((2*((abs(u[i]))^q)*((abs(u1[i]))^q) + tol)/((abs(u[i]))^(2*q) + (abs(u1[i]))^(2*q) + tol)));
        //y[i] = (1-tol)*(1 + tol - ((2*(noevent((abs(u[i]))^q))*(noevent((abs(u1[i]))^q)) + tol)/((noevent((abs(u[i]))^(2*q)) + noevent((abs(u1[i]))^(2*q)) + tol))));
        //y[i] = (1-tol)*(1 + tol - (2*(u[i]^(1.4))*(u1[i]^(1.4)) + tol)/(u[i]^(2.8) + u1[i]^(2.8) + tol));
        end for;

        // equation

        // for i in 1:n+1 loop
        //   y[i] = (1-tol)*(1 + tol - ((2*((abs(u[i]))^q)*((abs(u1[i]))^q) + tol)/((abs(u[i]))^(2*q) + (abs(u1[i]))^(2*q) + tol)));
        //   //y[i] = (1-tol)*(1 + tol - ((2*(noevent((abs(u[i]))^q))*(noevent((abs(u1[i]))^q)) + tol)/((noevent((abs(u[i]))^(2*q)) + noevent((abs(u1[i]))^(2*q)) + tol))));
        //   //y[i] = (1-tol)*(1 + tol - (2*(u[i]^(1.4))*(u1[i]^(1.4)) + tol)/(u[i]^(2.8) + u1[i]^(2.8) + tol));
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
            annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(
              Text(
                extent=[-34,40; 44,-30],
                style(color=3, rgbcolor={0,0,255}),
                string="c1"),
              Text(
                extent=[-100,72; -56,50],
                style(color=3, rgbcolor={0,0,255}),
                string="d1"),
              Text(
                extent=[-98,-48; -58,-70],
                style(color=3, rgbcolor={0,0,255}),
                string="d2")),
            Documentation(info="<html>
<p>
Implements the constant value c<sub>1</sub> needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log9.png\">

<p>
where <i>tol = 0.1h^q</i>, <i>q = 1.4</i> typically.
</p>

</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end C1;

        block C2
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        // parameter Integer gcl = worldModel1.gcl;
        // parameter Integer gcr = worldModel1.gcr;

        equation
        for i in 1:n+1 loop
          y[i] = u[i]/(u[i]-1);
        end for;

        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;
        // for i in gcl+1:gcl+n loop
        //   y[i] = u[i]/(u[i]-1);
        // end for;
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-42,38; 34,-30],
                style(color=3, rgbcolor={0,0,255}),
                string="c2"), Text(
                extent=[-92,12; -70,-10],
                style(color=3, rgbcolor={0,0,255}),
                string="c1")),
            Documentation(info="<html>
<p>
Implements the constant value c<sub>2</sub> needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log10.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end C2;

        block C3
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        // parameter Integer gcl = worldModel1.gcl;
        // parameter Integer gcr = worldModel1.gcr;

        equation
        for i in 1:n+1 loop
          y[i] = ((u2[i] - 1)*(u1[i]*(1 - u3[i]) - u[i]))/(u3[i] - u2[i]);
        end for;

        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;
        // for i in gcl+1:gcl+n loop
        //   y[i] = (u2[i] - 1)*(u1[i]*(1 - u3[i]) - u[i])/(u3[i] - u2[i]);
        // end for;
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
            annotation (extent=[-118,68; -100,92]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
            annotation (extent=[-118,28; -100,52]);
          Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
            annotation (extent=[-118,-52; -100,-28]);
          Modelica.Blocks.Interfaces.RealInput u3[worldModel1.n + 1]
            annotation (extent=[-118,-92; -100,-68]);
          annotation (Diagram, Icon(
              Text(
                extent=[-100,86; -68,72],
                style(color=3, rgbcolor={0,0,255}),
                string="d1"),
              Text(
                extent=[-100,46; -68,32],
                style(color=3, rgbcolor={0,0,255}),
                string="d2"),
              Text(
                extent=[-100,-34; -68,-48],
                style(color=3, rgbcolor={0,0,255}),
                string="c1"),
              Text(
                extent=[-100,-74; -66,-88],
                style(color=3, rgbcolor={0,0,255}),
                string="c2"),
              Text(
                extent=[-40,34; 40,-28],
                style(color=3, rgbcolor={0,0,255}),
                string="c3")),
            Documentation(info="<html>
<p>
Implements the constant value c<sub>3</sub> needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log11.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
        end C3;

        block C4
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        // parameter Integer gcl = worldModel1.gcl;
        // parameter Integer gcr = worldModel1.gcr;

        equation
        for i in 1:n+1 loop
          y[i] = u[i] - u1[i];
        end for;

        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;
        // for i in gcl+1:gcl+n loop
        //   y[i] = u[i] - u1[i];
        // end for;
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
            annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(
              Text(
                extent=[-100,68; -60,52],
                style(color=3, rgbcolor={0,0,255}),
                string="d1"),
              Text(
                extent=[-104,-52; -56,-68],
                style(color=3, rgbcolor={0,0,255}),
                string="c3"),
              Text(
                extent=[-26,38; 40,-30],
                style(color=3, rgbcolor={0,0,255}),
                string="c4")),
            Documentation(info="<html>
<p>
Implements the constant value c<sub>4</sub> needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log12.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end C4;
        annotation (Documentation(info="<html>
<p>
Implements the constant values c<sub>1</sub>, c<sub>2</sub>, c<sub>3</sub> and c<sub>4</sub> needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log3.png\">

<p>
where <i>tol = 0.1h^q</i>, <i>q = 1.4</i> typically.
</p>

</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end c;

      package n
        block n_plus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        // parameter Integer gcl = worldModel1.gcl;
        // parameter Integer gcr = worldModel1.gcr;

        equation
        for i in 1:n+1 loop
          y[i] = (-(log(1 - u[i])) - u[i])/((u[i])^2);
        end for;

        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;
        // for i in gcl+1:gcl+n loop
        //   y[i] = (-log(1 - u[i]) - u[i])/(u[i])^2;
        // end for;
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-22,36; 36,-24],
                style(color=3, rgbcolor={0,0,255}),
                string="n+")),
            Documentation(info="<html>
<p>
Implements the function n<sup>+</sup>(x) needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log13.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end n_plus;

        block n_minus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        // parameter Integer gcl = worldModel1.gcl;
        // parameter Integer gcr = worldModel1.gcr;

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
        equation

        for i in 1:n+1 loop
          y[i] = ((u[i] - 1)*(log(1 - u[i])) - u[i])/((u[i])^2);
        end for;

        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;
        // for i in gcl+1:gcl+n loop
        //   y[i] = -((u[i] - 1)*log(1 - u[i]) - u[i])/(u[i])^2;
        // end for;
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

          annotation (Diagram, Icon(Text(
                extent=[-32,36; 44,-26],
                style(color=3, rgbcolor={0,0,255}),
                string="n-")),
            Documentation(info="<html>
<p>
Implements the function n<sup>-</sup>(x) needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log14.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        end n_minus;
        annotation (Documentation(info="<html>
<p>
Implements the functions n<sup>+</sup>(x) and n<sup>-</sup>(x) needed for the logarithmic reconstruction
</p>

<img align=middle src=\"..\\Images\\log2.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end n;

      package R
        block Rminus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        inner parameter Integer n = worldModel1.n;
        inner parameter Integer gcl = worldModel1.gcl;
        inner parameter Integer gcr = worldModel1.gcr;

        // outer PDE.World.worldModel worldModel1;
        // parameter Integer n = worldModel1.n;
        // parameter Real deltax = 1/(n-1);
        //
        // equation
        // for i in 2:n+1 loop
        //   y[i] = u[i] + u1[i]*u3[i]*deltax + u2[i]*u4[i]*deltax;
        // end for;

          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
            annotation (extent=[-116,70; -100,90]);
          Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
            annotation (extent=[-116,50; -100,70]);
          Modelica.Blocks.Interfaces.RealInput u3[worldModel1.n + 1]
            annotation (extent=[-116,-70; -100,-50]);
          Modelica.Blocks.Interfaces.RealInput u4[worldModel1.n + 1]
            annotation (extent=[-116,-90; -100,-70]);
          annotation (Icon(
              Text(
                extent=[-98,86; -76,74],
                style(color=3, rgbcolor={0,0,255}),
                string="c1"),
              Text(
                extent=[-100,66; -72,54],
                style(color=3, rgbcolor={0,0,255}),
                string="c2"),
              Text(
                extent=[-98,-54; -74,-66],
                style(color=3, rgbcolor={0,0,255}),
                string="c3"),
              Text(
                extent=[-100,-74; -72,-86],
                style(color=3, rgbcolor={0,0,255}),
                string="c4"),
              Text(
                extent=[-34,34; 30,-30],
                style(color=3, rgbcolor={0,0,255}),
                string="R-")),
                            Diagram,
            Documentation(info="<html>
<p>
Implements
</p>

<img align=middle src=\"..\\Images\\log16.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          Modelica.Blocks.Math.Product product[worldModel1.n + 1]
            annotation (extent=[-46,72; -36,82]);
          Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
            annotation (extent=[-46,52; -36,62]);
          Modelica.Blocks.Math.Add add[worldModel1.n + 1]
            annotation (extent=[-20,62; -10,72]);
          Modelica.Blocks.Math.Product product2[worldModel1.n + 1]
            annotation (extent=[50,-6; 62,6]);
          Modelica.Blocks.Sources.RealExpression deltax[worldModel1.n + 1](y=1/(
                worldModel1.n))     annotation (extent=[-34,-26; -6,-6]);
          PDE.FiniteVolume.LDLR.n.n_plus n_plus1
                         annotation (extent=[-80,74; -68,86]);
          PDE.FiniteVolume.LDLR.n.n_plus n_plus2
                         annotation (extent=[-80,54; -68,66]);
        equation
          connect(u1, n_plus1.u) annotation (points=[-108,80; -81.2,80], style(
                color=74, rgbcolor={0,0,127}));
          connect(u2, n_plus2.u) annotation (points=[-108,60; -81.2,60], style(
                color=74, rgbcolor={0,0,127}));
          connect(n_plus1.y, product.u1) annotation (points=[-67.4,80; -47,80],
              style(color=74, rgbcolor={0,0,127}));
          connect(u3, product.u2) annotation (points=[-108,-60; -60,-60; -60,74;
                -47,74], style(color=74, rgbcolor={0,0,127}));
          connect(n_plus2.y, product1.u1) annotation (points=[-67.4,60; -47,60],
              style(color=74, rgbcolor={0,0,127}));
          connect(u4, product1.u2) annotation (points=[-108,-80; -54,-80; -54,54;
                -47,54], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, add.u1) annotation (points=[-35.5,77; -26,77; -26,70;
                -21,70], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, add.u2) annotation (points=[-35.5,57; -26,57; -26,
                64; -21,64], style(color=74, rgbcolor={0,0,127}));
          connect(deltax.y, product2.u2) annotation (points=[-4.6,-16; 40,-16;
                40,-3.6; 48.8,-3.6],
                                  style(color=74, rgbcolor={0,0,127}));
          connect(add.y, product2.u1) annotation (points=[-9.5,67; 40,67; 40,
                3.6; 48.8,3.6],
                           style(color=74, rgbcolor={0,0,127}));
          connect(product2.y, y) annotation (points=[62.6,0; 110,0], style(color=
                  74, rgbcolor={0,0,127}));
        end Rminus;

        block Rplus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        inner parameter Integer n = worldModel1.n;
        inner parameter Integer gcl = worldModel1.gcl;
        inner parameter Integer gcr = worldModel1.gcr;

          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
            annotation (extent=[-116,68; -100,92]);
          Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + 1]
            annotation (extent=[-116,38; -100,62]);
          Modelica.Blocks.Interfaces.RealInput u3[worldModel1.n + 1]
            annotation (extent=[-116,-62; -100,-38]);
          Modelica.Blocks.Interfaces.RealInput u4[worldModel1.n + 1]
            annotation (extent=[-116,-92; -100,-68]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(
              Text(
                extent=[-100,86; -70,72],
                style(color=3, rgbcolor={0,0,255}),
                string="c1"),
              Text(
                extent=[-96,58; -72,44],
                style(color=3, rgbcolor={0,0,255}),
                string="c2"),
              Text(
                extent=[-102,-44; -66,-58],
                style(color=3, rgbcolor={0,0,255}),
                string="c3"),
              Text(
                extent=[-104,-72; -64,-86],
                style(color=3, rgbcolor={0,0,255}),
                string="c4"),
              Text(
                extent=[-30,34; 34,-32],
                style(color=3, rgbcolor={0,0,255}),
                string="R+")),
            Documentation(info="<html>
<p>
Implements
</p>

<img align=middle src=\"..\\Images\\log17.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          PDE.FiniteVolume.LDLR.n.n_minus n_minus1
                           annotation (extent=[-80,74; -68,86]);
          PDE.FiniteVolume.LDLR.n.n_minus n_minus2
                           annotation (extent=[-80,44; -68,56]);
          Modelica.Blocks.Math.Product product[worldModel1.n + 1]
            annotation (extent=[-50,74; -38,86]);
          Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
            annotation (extent=[-50,44; -38,56]);
          Modelica.Blocks.Math.Add add[worldModel1.n + 1]
            annotation (extent=[-20,58; -8,70]);
          Modelica.Blocks.Math.Product product2[worldModel1.n + 1]
            annotation (extent=[48,-6; 60,6]);
          Modelica.Blocks.Sources.RealExpression deltax[worldModel1.n + 1](y=1/(
                worldModel1.n))     annotation (extent=[-34,-28; -6,-8]);
        equation
          connect(u1, n_minus1.u) annotation (points=[-108,80; -81.2,80], style(
                color=74, rgbcolor={0,0,127}));
          connect(u2, n_minus2.u) annotation (points=[-108,50; -81.2,50], style(
                color=74, rgbcolor={0,0,127}));
          connect(n_minus1.y, product.u1) annotation (points=[-67.4,80; -60,80;
                -60,83.6; -51.2,83.6], style(color=74, rgbcolor={0,0,127}));
          connect(n_minus2.y, product1.u1) annotation (points=[-67.4,50; -60,50;
                -60,53.6; -51.2,53.6], style(color=74, rgbcolor={0,0,127}));
          connect(u3, product.u2) annotation (points=[-108,-50; -58,-50; -58,76.4;
                -51.2,76.4], style(color=74, rgbcolor={0,0,127}));
          connect(u4, product1.u2) annotation (points=[-108,-80; -56,-80; -56,
                46.4; -51.2,46.4], style(color=74, rgbcolor={0,0,127}));
          connect(product.y, add.u1) annotation (points=[-37.4,80; -30,80; -30,
                67.6; -21.2,67.6], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, add.u2) annotation (points=[-37.4,50; -30,50; -30,
                60.4; -21.2,60.4], style(color=74, rgbcolor={0,0,127}));
          connect(add.y, product2.u1) annotation (points=[-7.4,64; 20,64; 20,
                3.6; 46.8,3.6],
                           style(color=74, rgbcolor={0,0,127}));
          connect(deltax.y, product2.u2) annotation (points=[-4.6,-18; 20,-18;
                20,-4; 46.8,-4; 46.8,-3.6],
                                         style(color=74, rgbcolor={0,0,127}));
          connect(product2.y, y) annotation (points=[60.6,0; 110,0], style(color=
                  74, rgbcolor={0,0,127}));
        end Rplus;
        annotation (Documentation(info="<html>
<p>
Implements
</p>

<img align=middle src=\"..\\Images\\log15.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end R;

      package L
        block LDLRminus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        inner parameter Integer n = worldModel1.n;
        inner parameter Integer gcl = worldModel1.gcl;
        inner parameter Integer gcr = worldModel1.gcr;

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-40,28; 64,-26],
                style(color=3, rgbcolor={0,0,255}),
                string="LDLR - ")),
            Documentation(info="<html>
<p>
Combine d<sub>1</sub>-, d<sub>2</sub>-, c<sub>1</sub>, c<sub>2</sub>, c<sub>3</sub>, c<sub>4</sub> and R<sup>-</sup> blocks to implement
</p>

<img align=middle src=\"..\\Images\\log16.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          PDE.FiniteVolume.LDLR.d.D1minus d1_1
                                             annotation (extent=[-80,14; -68,26]);
          PDE.FiniteVolume.LDLR.d.D2minus d2_1
            annotation (extent=[-80,-26; -68,-14]);
          PDE.FiniteVolume.LDLR.c.C1 c1_1
                  annotation (extent=[-48,-26; -34,-14]);
          PDE.FiniteVolume.LDLR.c.C2 c2_1
                  annotation (extent=[-22,10; -10,22]);
          PDE.FiniteVolume.LDLR.c.C3 c3_1
                  annotation (extent=[0,-16; 18,-4]);
          PDE.FiniteVolume.LDLR.c.C4 c4_1
                  annotation (extent=[28,-26; 46,-14]);
          PDE.FiniteVolume.LDLR.R.Rminus reconstruction
            annotation (extent=[56,-14; 88,14]);
        equation
          connect(u, d1_1.u) annotation (points=[-120,0; -90,0; -90,20; -81.2,20],
              style(color=74, rgbcolor={0,0,127}));
          connect(u, d2_1.u) annotation (points=[-120,0; -90,0; -90,-20; -81.2,
                -20], style(color=74, rgbcolor={0,0,127}));
          connect(d1_1.y, c1_1.u) annotation (points=[-67.4,20; -60,20; -60,-16.4;
                -49.4,-16.4], style(color=74, rgbcolor={0,0,127}));
          connect(d2_1.y, c1_1.u1) annotation (points=[-67.4,-20; -60,-20; -60,
                -23.6; -49.4,-23.6], style(color=74, rgbcolor={0,0,127}));
          connect(c1_1.y, c2_1.u) annotation (points=[-33.3,-20; -28,-20; -28,16;
                -23.2,16], style(color=74, rgbcolor={0,0,127}));
          connect(c1_1.y, c3_1.u2) annotation (points=[-33.3,-20; -6,-20; -6,
                -12.4; -0.81,-12.4], style(color=74, rgbcolor={0,0,127}));
          connect(c2_1.y, c3_1.u3) annotation (points=[-9.4,16; -4,16; -4,-14.8;
                -0.81,-14.8], style(color=74, rgbcolor={0,0,127}));
          connect(d1_1.y, c3_1.u) annotation (points=[-67.4,20; -48,20; -48,
                -5.2; -0.81,-5.2],
                             style(color=74, rgbcolor={0,0,127}));
          connect(d2_1.y, c3_1.u1) annotation (points=[-67.4,-20; -64,-20; -64,
                -7.6; -0.81,-7.6], style(color=74, rgbcolor={0,0,127}));
          connect(d1_1.y, c4_1.u) annotation (points=[-67.4,20; -40,20; -40,4;
                22,4; 22,-16.4; 26.2,-16.4],
                                          style(color=74, rgbcolor={0,0,127}));
          connect(c3_1.y, c4_1.u1) annotation (points=[18.9,-10; 20,-10; 20,
                -23.6; 26.2,-23.6],
                             style(color=74, rgbcolor={0,0,127}));
          connect(c4_1.y, reconstruction.u4) annotation (points=[46.9,-20; 50,-20;
                50,-11.2; 54.72,-11.2], style(color=74, rgbcolor={0,0,127}));
          connect(c3_1.y, reconstruction.u3) annotation (points=[18.9,-10; 46,
                -10; 46,-8.4; 54.72,-8.4],
                                      style(color=74, rgbcolor={0,0,127}));
          connect(c2_1.y, reconstruction.u2) annotation (points=[-9.4,16; 40,16;
                40,8.4; 54.72,8.4], style(color=74, rgbcolor={0,0,127}));
          connect(c1_1.y, reconstruction.u1) annotation (points=[-33.3,-20; -28,
                -20; -28,24; 46,24; 46,11.2; 54.72,11.2], style(color=74,
                rgbcolor={0,0,127}));
          connect(reconstruction.y, y) annotation (points=[89.6,1.77636e-015;
                96.8,1.77636e-015; 96.8,0; 110,0], style(color=74, rgbcolor={0,0,
                  127}));
        end LDLRminus;

        block LDLRplus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        inner parameter Integer n = worldModel1.n;
        inner parameter Integer gcl = worldModel1.gcl;
        inner parameter Integer gcr = worldModel1.gcr;

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-44,26; 66,-24],
                style(color=3, rgbcolor={0,0,255}),
                string="LDLR + ")),
            Documentation(info="<html>
<p>
Combine d<sub>1</sub>+, d<sub>2</sub>+, c<sub>1</sub>, c<sub>2</sub>, c<sub>3</sub>, c<sub>4</sub> and R<sup>+</sup> blocks to implement
</p>

<img align=middle src=\"..\\Images\\log17.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          PDE.FiniteVolume.LDLR.d.D1plus d1plus
                        annotation (extent=[-80,12; -64,28]);
          PDE.FiniteVolume.LDLR.d.D2plus d2plus
                        annotation (extent=[-80,-28; -64,-12]);
          PDE.FiniteVolume.LDLR.c.C1 c1_1
                  annotation (extent=[-40,-26; -26,-14]);
          PDE.FiniteVolume.LDLR.c.C2 c2_1
                  annotation (extent=[-14,14; -2,26]);
          PDE.FiniteVolume.LDLR.c.C3 c3_1
                  annotation (extent=[6,-12; 24,2]);
          PDE.FiniteVolume.LDLR.c.C4 c4_1
                  annotation (extent=[36,-30; 50,-18]);
          PDE.FiniteVolume.LDLR.R.Rplus reconstructionPlus
            annotation (extent=[58,-12; 84,12]);
        equation
          connect(u, d1plus.u) annotation (points=[-120,0; -92,0; -92,20; -81.6,
                20], style(color=74, rgbcolor={0,0,127}));
          connect(u, d2plus.u) annotation (points=[-120,0; -92,0; -92,-20; -81.6,
                -20], style(color=74, rgbcolor={0,0,127}));
          connect(d1plus.y, c1_1.u) annotation (points=[-63.2,20; -50,20; -50,
                -16.4; -41.4,-16.4], style(color=74, rgbcolor={0,0,127}));
          connect(d2plus.y, c1_1.u1) annotation (points=[-63.2,-20; -52,-20; -52,
                -23.6; -41.4,-23.6], style(color=74, rgbcolor={0,0,127}));
          connect(c1_1.y, c2_1.u) annotation (points=[-25.3,-20; -20,-20; -20,
                20; -15.2,20],
                           style(color=74, rgbcolor={0,0,127}));
          connect(d1plus.y, c3_1.u) annotation (points=[-63.2,20; -40,20; -40,
                0.6; 5.19,0.6],
                           style(color=74, rgbcolor={0,0,127}));
          connect(d2plus.y, c3_1.u1) annotation (points=[-63.2,-20; -52,-20;
                -52,-2.2; 5.19,-2.2],
                                  style(color=74, rgbcolor={0,0,127}));
          connect(c1_1.y, c3_1.u2) annotation (points=[-25.3,-20; -6,-20; -6,
                -7.8; 5.19,-7.8],
                            style(color=74, rgbcolor={0,0,127}));
          connect(c2_1.y, c3_1.u3) annotation (points=[-1.4,20; 0,20; 0,-10.6;
                5.19,-10.6], style(color=74, rgbcolor={0,0,127}));
          connect(d1plus.y, c4_1.u) annotation (points=[-63.2,20; -32,20; -32,6;
                -12,6; -12,-18; 20,-18; 20,-20.4; 34.6,-20.4], style(color=74,
                rgbcolor={0,0,127}));
          connect(c3_1.y, c4_1.u1) annotation (points=[24.9,-5; 28,-5; 28,-27.6;
                34.6,-27.6], style(color=74, rgbcolor={0,0,127}));
          connect(c4_1.y, reconstructionPlus.u4) annotation (points=[50.7,-24; 54,
                -24; 54,-9.6; 56.96,-9.6], style(color=74, rgbcolor={0,0,127}));
          connect(c3_1.y, reconstructionPlus.u3) annotation (points=[24.9,-5;
                40.45,-5; 40.45,-6; 56.96,-6], style(color=74, rgbcolor={0,0,127}));
          connect(c2_1.y, reconstructionPlus.u2) annotation (points=[-1.4,20;
                44,20; 44,6; 56.96,6],
                                    style(color=74, rgbcolor={0,0,127}));
          connect(c1_1.y, reconstructionPlus.u1) annotation (points=[-25.3,-20;
                -22,-20; -22,32; 50,32; 50,9.6; 56.96,9.6], style(color=74,
                rgbcolor={0,0,127}));
          connect(reconstructionPlus.y, y) annotation (points=[85.3,0; 110,0],
              style(color=74, rgbcolor={0,0,127}));
        end LDLRplus;

        annotation (Documentation(info="<html>
<p>
This package contains two blocks, LDLRminus and LDLRplus, that combine d<sub>1</sub>-, d<sub>2</sub>-, d<sub>1</sub>+, d<sub>2</sub>+, c<sub>1</sub>, c<sub>2</sub>, c<sub>3</sub>, c<sub>4</sub> and R<sup>-</sup> blocks to implement
</p>

<p>
<img align=middle src=\"..\\Images\\log16.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\log17.png\">
</p>


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end L;

      package u
        block u_minus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;

        equation
        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;
        for i in 1:n+1 loop
          y[i] = u[i+gcl-1] + u1[i];
        end for;
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

          annotation (Icon(Text(
                extent=[-24,28; 56,-22],
                style(color=3, rgbcolor={0,0,255}),
                string="u - ")),  Diagram,
            Documentation(info="<html>
<p>
Implements the reconstruction of the value u<sup>-</sup> at the interface x<sub>i+1/2</sub>
</p>

<img align=middle src=\"..\\Images\\log18.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
        end u_minus;

        block u_plus
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;

        equation
        for i in 1:n+1 loop
          y[i] = u[i+gcl] + u1[i];
        end for;

        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;
        // for i in gcl+1:gcl+n loop
        //   y[i] = u[i+1] + u1[i];
        // end for;
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

          annotation (Icon(Text(
                extent=[-22,30; 56,-20],
                style(color=3, rgbcolor={0,0,255}),
                string="u + ")),   Diagram,
            Documentation(info="<html>
<p>
Implements the reconstruction of the value u<sup>+</sup> at the interface x<sub>i+1/2</sub>
</p>

<img align=middle src=\"..\\Images\\log19.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl
             + worldModel1.gcr]
            annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
        end u_plus;

        annotation (Documentation(info="<html>
<p>
Implements the reconstruction of the values u<sup>-</sup> and u<sup>+</sup> at the interface x<sub>i+1/2</sub>
</p>

<img align=middle src=\"..\\Images\\log1.png\">


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end u;

      annotation (Documentation(info="<html>
<p>
Implements the limiter-free third order logarithmic reconstruction (LDLR). In LDLR the values of the flux at the interface x<sub>i+1/2</sub> are <br>
reconstructed in the following way
</p>

<img align=middle src=\"..\\Images\\log1.png\">

<p>
where
</p>

<p>
<img align=middle src=\"..\\Images\\log2.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\log3.png\">
</p>

<p>
where <i>tol = 0.1h^q</i>, <i>q = 1.4</i> typically and the lateral derivatives d<sub>1</sub> and d<sub>2</sub> are obtained in the following way
</p>

<img align=middle src=\"..\\Images\\log4.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      package Reconstruction
        block Rec
          import PDE;
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
          parameter Integer method = worldModel1.qss;
          inner parameter Integer n = worldModel1.n;
          inner parameter Integer gcl = worldModel1.gcl;
          inner parameter Integer gcr = worldModel1.gcr;

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.
            gcl + worldModel1.gcr] annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
            annotation (extent=[100,50; 120,70]);
          Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.n + 1]
            annotation (extent=[100,-70; 120,-50]);
          annotation (Diagram, Icon(
              Text(
                extent=[-28,32; 22,-24],
                style(color=3, rgbcolor={0,0,255}),
                string="u"),
              Text(
                extent=[66,74; 92,48],
                style(color=3, rgbcolor={0,0,255}),
                string="-"),
              Text(
                extent=[64,-50; 94,-70],
                style(color=3, rgbcolor={0,0,255}),
                string="+"),
              Text(
                extent=[-94,10; -68,-12],
                style(color=3, rgbcolor={0,0,255}),
                string="Q")),
            Documentation(info="<html>
<p>
Implements the reconstruction of the values u<sup>-</sup> and u<sup>+</sup> at the interface x<sub>i+1/2</sub> by using the
LDLR-, LDLR+, u- and u+ blocks.
</p>

<p>
<img align=middle src=\"..\\Images\\log18.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\log19.png\">
</p>


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
          PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus
                                          annotation (extent=[-40,20; -20,40]);
          PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus
                                        annotation (extent=[-40,-40; -20,-20]);
          PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus
            annotation (extent=[0,20; 20,40]);
          PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus
            annotation (extent=[0,-40; 20,-20]);
        equation
          connect(lDLRminus_plus.y, u_minus_plus.u1) annotation (points=[-19,30;
                -12,30; -12,24; -2,24], style(color=74, rgbcolor={0,0,127}));
          connect(lDLRplus_plus.y, u_plus_plus.u1) annotation (points=[-19,-30;
                -12,-30; -12,-36; -2,-36], style(color=74, rgbcolor={0,0,127}));
          connect(u, lDLRminus_plus.u) annotation (points=[-120,0; -60,0; -60,
                30; -42,30], style(color=74, rgbcolor={0,0,127}));
          connect(u, lDLRplus_plus.u) annotation (points=[-120,0; -60,0; -60,
                -30; -42,-30], style(color=74, rgbcolor={0,0,127}));
          connect(u, u_minus_plus.u) annotation (points=[-120,0; -60,0; -60,46;
                -12,46; -12,36; -2,36], style(color=74, rgbcolor={0,0,127}));
          connect(u, u_plus_plus.u) annotation (points=[-120,0; -60,0; -60,-14;
                -12,-14; -12,-24; -2,-24], style(color=74, rgbcolor={0,0,127}));
          connect(u_minus_plus.y, y) annotation (points=[21,30; 70,30; 70,60;
                110,60], style(color=74, rgbcolor={0,0,127}));
          connect(u_plus_plus.y, y1) annotation (points=[21,-30; 70,-30; 70,-60;
                110,-60], style(color=74, rgbcolor={0,0,127}));
        end Rec;
        annotation (Documentation(info="<html>
<p>
This package contains block Rec that computes the the values u<sup>-</sup> and u<sup>+</sup> at the interface x<sub>i+1/2</sub> by using the
LDLR-, LDLR+, u- and u+ blocks.
</p>

<p>
<img align=middle src=\"..\\Images\\log18.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\log19.png\">
</p>


</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end Reconstruction;
    end LDLR;

    package Examples

      package Burger

        model BurgerEquationLDLR
          PDE.FiniteVolume.FVMIntegrator.FVIntegrator Burger(
                                                      bcr=0,
            gcl=2,
            gcr=2,
            vb=worldModel1.gcl + 1,
            icb=worldModel1.gcl + 1,
            bcl=0)
            annotation (extent=[-60,24; -20,50]);
          PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus
                                             annotation (extent=[0,52; 22,60]);
          PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus
                                           annotation (extent=[0,40; 22,48]);
          PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus
                                         annotation (extent=[40,52; 56,60]);
          PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus
                                       annotation (extent=[40,40; 56,48]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Burger큦 equation</font></h3>
<p>
Implements the inviscid Burger큦 equation
</p>

<img align=middle src=\"..\\Images\\b1.png\">

<p>
The initial condition is
</p>

<img align=middle src=\"..\\Images\\b2.png\">

<p>
and boundary conditions are
</p>

<p>
<img align=middle src=\"..\\Images\\b3.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\b4.png\">
</p>

<p>
The analytical solution of this problem is implemented in PDE->MOL->Examples->Burger->BAN block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Product product[worldModel1.n + 1]
                             annotation (extent=[-50,-4; -44,2]);
          Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
                                annotation (extent=[-50,-16; -44,-10]);
          Modelica.Blocks.Math.Division division[worldModel1.n + 1]
                                annotation (extent=[-26,-4; -20,2]);
          Modelica.Blocks.Math.Division division1[worldModel1.n + 1]
                                annotation (extent=[-26,-16; -20,-10]);
          Modelica.Blocks.Sources.RealExpression const[worldModel1.n + 1](y=2.0)
            annotation (extent=[-50,-34; -44,-20]);
          PDE.FiniteVolume.Fluxes.LaxFriedrichFlux.LF lF(
                                              alpha=0.9)
                                annotation (extent=[0,-16; 18,2]);
          inner World.worldModel worldModel1(
            gcl=2,
            gcr=2,
            n=40)
            annotation (extent=[-20,88; 20,100]);
          MOL.Examples.Burger.BICx bICx annotation (extent=[-88,32; -78,42]);
          Modelica.Blocks.Sources.RealExpression BC[2]
            annotation (extent=[-92,10; -74,30]);
        equation
          connect(Burger.y, lDLRminus_plus.u)    annotation (points=[-18,44.8;
                -12,44.8; -12,56; -2.2,56], style(color=74, rgbcolor={0,0,127}));
          connect(Burger.y, lDLRplus_plus.u)    annotation (points=[-18,44.8;
                -12,44.8; -12,44; -2.2,44],
                                        style(color=74, rgbcolor={0,0,127}));
          connect(Burger.y, u_minus_plus.u)    annotation (points=[-18,44.8; -12,
                44.8; -12,66; 32,66; 32,58.4; 38.4,58.4], style(color=74,
                rgbcolor={0,0,127}));
          connect(Burger.y, u_plus_plus.u)    annotation (points=[-18,44.8; -12,
                44.8; -12,50; 32,50; 32,46.4; 38.4,46.4], style(color=74,
                rgbcolor={0,0,127}));
          connect(lDLRminus_plus.y, u_minus_plus.u1) annotation (points=[23.1,56;
                32,56; 32,53.6; 38.4,53.6], style(color=74, rgbcolor={0,0,127}));
          connect(lDLRplus_plus.y, u_plus_plus.u1) annotation (points=[23.1,44;
                32,44; 32,41.6; 38.4,41.6], style(color=74, rgbcolor={0,0,127}));
          connect(u_minus_plus.y, product.u1) annotation (points=[56.8,56; 80,
                56; 80,6; -54,6; -54,0.8; -50.6,0.8],
                                            style(color=74, rgbcolor={0,0,127}));
          connect(u_minus_plus.y, product.u2) annotation (points=[56.8,56; 80,
                56; 80,6; -54,6; -54,-2.8; -50.6,-2.8],
                                                    style(color=74, rgbcolor={0,0,
                  127}));
          connect(u_plus_plus.y, product1.u1) annotation (points=[56.8,44; 78,
                44; 78,8; -56,8; -56,-11.2; -50.6,-11.2],
                                                      style(color=74, rgbcolor={0,
                  0,127}));
          connect(u_plus_plus.y, product1.u2) annotation (points=[56.8,44; 78,
                44; 78,8; -56,8; -56,-14.8; -50.6,-14.8],
                                                      style(color=74, rgbcolor={0,
                  0,127}));
          connect(product.y, division.u1) annotation (points=[-43.7,-1; -36,-1;
                -36,0.8; -26.6,0.8], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, division1.u1) annotation (points=[-43.7,-13; -36,
                -13; -36,-11.2; -26.6,-11.2], style(color=74, rgbcolor={0,0,127}));
          connect(const.y, division.u2) annotation (points=[-43.7,-27; -32,-27;
                -32,-2.8; -26.6,-2.8], style(color=74, rgbcolor={0,0,127}));
          connect(const.y, division1.u2) annotation (points=[-43.7,-27; -32,-27;
                -32,-14.8; -26.6,-14.8], style(color=74, rgbcolor={0,0,127}));
          connect(division.y, lF.u) annotation (points=[-19.7,-1; -8,-1; -8,0.2;
                -0.81,0.2], style(color=74, rgbcolor={0,0,127}));
          connect(division1.y, lF.u1) annotation (points=[-19.7,-13; -8,-13; -8,
                -2.5; -0.81,-2.5], style(color=74, rgbcolor={0,0,127}));
          connect(u_plus_plus.y, lF.u2) annotation (points=[56.8,44; 78,44; 78,
                8; -6,8; -6,-11.5; -0.81,-11.5],
                                              style(color=74, rgbcolor={0,0,127}));
          connect(u_minus_plus.y, lF.u3) annotation (points=[56.8,56; 80,56; 80,
                6; -4,6; -4,-14.2; -0.81,-14.2],
                                              style(color=74, rgbcolor={0,0,127}));
          connect(lF.y, Burger.u)          annotation (points=[18.9,-7; 88,-7; 88,
                72; -70,72; -70,47.4; -62.2,47.4], style(color=74, rgbcolor={0,0,
                  127}));
          connect(BC.y, Burger.u3) annotation (points=[-73.1,20; -68,20; -68,
                30.5; -62.2,30.5], style(color=74, rgbcolor={0,0,127}));
          connect(BC.y, Burger.u4) annotation (points=[-73.1,20; -68,20; -68,
                26.6; -62.2,26.6], style(color=74, rgbcolor={0,0,127}));
          connect(bICx.y, Burger.u2) annotation (points=[-77.5,37; -70.75,37;
                -70.75,37; -62.2,37], style(color=74, rgbcolor={0,0,127}));
        end BurgerEquationLDLR;
        annotation (Documentation(info="<html>
<p>
This package contains Burger큦 equation solved with the Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Burger;

      package BuckleyLeverett
        model BuckleyLeverettEquation
          PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus
                                             annotation (extent=[-8,52; 14,60]);
          PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus
                                           annotation (extent=[-8,20; 14,28]);
          PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus
                                         annotation (extent=[26,52; 42,60]);
          PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus
                                       annotation (extent=[26,20; 42,28]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Buckley-Leverett equation</font></h3>
<p>
Implements the Buckley-Leverett equation
</p>

<img align=middle src=\"..\\Images\\bl1.png\">

<p>
The initial condition is
</p>

<img align=middle src=\"..\\Images\\bl2.png\">

<p>
and boundary conditions are
</p>

<p>
<img align=middle src=\"..\\Images\\bl3.png\">
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Math.Product product[worldModel1.n + 1]
            annotation (extent=[-64,-20; -56,-12]);
          Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
            annotation (extent=[-42,-20; -34,-12]);
          Modelica.Blocks.Math.Product product2[worldModel1.n + 1]
            annotation (extent=[-64,-50; -56,-42]);
          Modelica.Blocks.Math.Product product3[worldModel1.n + 1]
            annotation (extent=[-42,-50; -34,-42]);
          Modelica.Blocks.Sources.RealExpression four[worldModel1.n + 1](y=4)
            annotation (extent=[-64,-38; -56,-22]);
          Modelica.Blocks.Sources.RealExpression one[worldModel1.n + 1](y=1)
            annotation (extent=[-22,-38; -14,-22]);
          Modelica.Blocks.Math.Add add[worldModel1.n + 1](k2=-1)
            annotation (extent=[-4,-28; 4,-20]);
          Modelica.Blocks.Math.Add add1[worldModel1.n + 1](k2=-1)
            annotation (extent=[-4,-42; 4,-34]);
          Modelica.Blocks.Math.Product product4[worldModel1.n + 1]
            annotation (extent=[18,-28; 26,-20]);
          Modelica.Blocks.Math.Product product5[worldModel1.n + 1]
            annotation (extent=[18,-42; 26,-34]);
          Modelica.Blocks.Math.Add add2[worldModel1.n + 1]
            annotation (extent=[38,-28; 46,-20]);
          Modelica.Blocks.Math.Add add3[worldModel1.n + 1]
            annotation (extent=[38,-42; 46,-34]);
          Modelica.Blocks.Math.Division division[worldModel1.n + 1]
            annotation (extent=[56,-28; 64,-20]);
          Modelica.Blocks.Math.Division division1[worldModel1.n + 1]
            annotation (extent=[56,-42; 64,-34]);
          PDE.FiniteVolume.Fluxes.LaxFriedrichFlux.LF lF(
                                              alpha=0.5)
                     annotation (extent=[70,-36; 82,-24]);
          Modelica.Blocks.Sources.RealExpression IC[worldModel1.n](y=1)
            annotation (extent=[-88,34; -80,46]);
          Modelica.Blocks.Sources.RealExpression BCL[worldModel1.gcl]
            annotation (extent=[-88,24; -80,36]);
          Modelica.Blocks.Sources.RealExpression BCR[worldModel1.gcr]
            annotation (extent=[-88,14; -80,26]);
          inner World.worldModel worldModel1(n=40)
            annotation (extent=[-20,88; 20,100]);
          PDE.FiniteVolume.FVMIntegrator.FVIntegrator BuckleyLeverett(
            ve=worldModel1.gcl + worldModel1.n,
            ice=worldModel1.gcl + worldModel1.n,
            vb=worldModel1.gcl + 1,
            icb=worldModel1.gcl + 1,
            bcl=0,
            bcr=0)
            annotation (extent=[-66,20; -28,60]);
        equation
          connect(lDLRplus_plus.y, u_plus_plus.u1) annotation (points=[15.1,24;
                20,24; 20,21.6; 24.4,21.6], style(color=74, rgbcolor={0,0,127}));
          connect(lDLRminus_plus.y, u_minus_plus.u1) annotation (points=[15.1,56;
                20,56; 20,53.6; 24.4,53.6],     style(color=74, rgbcolor={0,0,
                  127}));
          connect(u_minus_plus.y, product.u1) annotation (points=[42.8,56; 66,
                56; 66,-6; -68,-6; -68,-13.6; -64.8,-13.6], style(color=74,
                rgbcolor={0,0,127}));
          connect(u_minus_plus.y, product.u2) annotation (points=[42.8,56; 66,
                56; 66,-6; -68,-6; -68,-18.4; -64.8,-18.4], style(color=74,
                rgbcolor={0,0,127}));
          connect(product.y, product1.u1) annotation (points=[-55.6,-16; -48,
                -16; -48,-13.6; -42.8,-13.6], style(color=74, rgbcolor={0,0,127}));
          connect(four.y, product1.u2) annotation (points=[-55.6,-30; -48,-30;
                -48,-18.4; -42.8,-18.4], style(color=74, rgbcolor={0,0,127}));
          connect(four.y, product3.u1) annotation (points=[-55.6,-30; -48,-30;
                -48,-43.6; -42.8,-43.6], style(color=74, rgbcolor={0,0,127}));
          connect(u_plus_plus.y, product2.u1) annotation (points=[42.8,24; 60,
                24; 60,0; -74,0; -74,-43.6; -64.8,-43.6], style(color=74,
                rgbcolor={0,0,127}));
          connect(u_plus_plus.y, product2.u2) annotation (points=[42.8,24; 60,
                24; 60,0; -74,0; -74,-48.4; -64.8,-48.4], style(color=74,
                rgbcolor={0,0,127}));
          connect(product2.y, product3.u2) annotation (points=[-55.6,-46; -48,
                -46; -48,-48.4; -42.8,-48.4], style(color=74, rgbcolor={0,0,127}));
          connect(one.y, add.u1) annotation (points=[-13.6,-30; -8,-30; -8,
                -21.6; -4.8,-21.6],
                            style(color=74, rgbcolor={0,0,127}));
          connect(one.y, add1.u1) annotation (points=[-13.6,-30; -8,-30; -8,
                -35.6; -4.8,-35.6],
                            style(color=74, rgbcolor={0,0,127}));
          connect(u_minus_plus.y, add.u2) annotation (points=[42.8,56; 66,56;
                66,-6; -10,-6; -10,-26.4; -4.8,-26.4],
                                                  style(color=74, rgbcolor={0,0,
                  127}));
          connect(u_plus_plus.y, add1.u2) annotation (points=[42.8,24; 46,24;
                46,0; -26,0; -26,-40.4; -4.8,-40.4],
                                                style(color=74, rgbcolor={0,0,
                  127}));
          connect(add.y, product4.u1) annotation (points=[4.4,-24; 12,-24; 12,
                -21.6; 17.2,-21.6], style(color=74, rgbcolor={0,0,127}));
          connect(add.y, product4.u2) annotation (points=[4.4,-24; 12,-24; 12,
                -26.4; 17.2,-26.4], style(color=74, rgbcolor={0,0,127}));
          connect(add1.y, product5.u1) annotation (points=[4.4,-38; 12,-38; 12,
                -35.6; 17.2,-35.6], style(color=74, rgbcolor={0,0,127}));
          connect(add1.y, product5.u2) annotation (points=[4.4,-38; 12,-38; 12,
                -40.4; 17.2,-40.4], style(color=74, rgbcolor={0,0,127}));
          connect(product4.y, add2.u2) annotation (points=[26.4,-24; 32,-24; 32,
                -26.4; 37.2,-26.4], style(color=74, rgbcolor={0,0,127}));
          connect(product5.y, add3.u1) annotation (points=[26.4,-38; 32,-38; 32,
                -35.6; 37.2,-35.6], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, add2.u1) annotation (points=[-33.6,-16; 32,-16;
                32,-21.6; 37.2,-21.6], style(color=74, rgbcolor={0,0,127}));
          connect(product3.y, add3.u2) annotation (points=[-33.6,-46; 32,-46;
                32,-40.4; 37.2,-40.4], style(color=74, rgbcolor={0,0,127}));
          connect(add2.y, division.u2) annotation (points=[46.4,-24; 50,-24; 50,
                -26.4; 55.2,-26.4], style(color=74, rgbcolor={0,0,127}));
          connect(add3.y, division1.u2) annotation (points=[46.4,-38; 50,-38;
                50,-40.4; 55.2,-40.4], style(color=74, rgbcolor={0,0,127}));
          connect(product1.y, division.u1) annotation (points=[-33.6,-16; 50,
                -16; 50,-21.6; 55.2,-21.6], style(color=74, rgbcolor={0,0,127}));
          connect(product3.y, division1.u1) annotation (points=[-33.6,-46; 52,
                -46; 52,-35.6; 55.2,-35.6], style(color=74, rgbcolor={0,0,127}));
          connect(division.y, lF.u) annotation (points=[64.4,-24; 66.93,-24;
                66.93,-25.2; 69.46,-25.2], style(color=74, rgbcolor={0,0,127}));
          connect(division1.y, lF.u1) annotation (points=[64.4,-38; 66,-38; 66,
                -27; 69.46,-27], style(color=74, rgbcolor={0,0,127}));
          connect(u_plus_plus.y, lF.u2) annotation (points=[42.8,24; 68,24; 68,
                -33; 69.46,-33], style(color=74, rgbcolor={0,0,127}));
          connect(u_minus_plus.y, lF.u3) annotation (points=[42.8,56; 66,56; 66,
                -34.8; 69.46,-34.8], style(color=74, rgbcolor={0,0,127}));
          connect(BuckleyLeverett.y, lDLRminus_plus.u) annotation (points=[-26.1,52;
                -18,52; -18,56; -10.2,56], style(color=74, rgbcolor={0,0,127}));
          connect(BuckleyLeverett.y, lDLRplus_plus.u) annotation (points=[-26.1,52;
                -18,52; -18,24; -10.2,24], style(color=74, rgbcolor={0,0,127}));
          connect(BuckleyLeverett.y, u_minus_plus.u) annotation (points=[-26.1,52;
                -18,52; -18,64; 20,64; 20,58.4; 24.4,58.4],   style(color=74,
                rgbcolor={0,0,127}));
          connect(BuckleyLeverett.y, u_plus_plus.u) annotation (points=[-26.1,52;
                -18,52; -18,32; 20,32; 20,26.4; 24.4,26.4],   style(color=74,
                rgbcolor={0,0,127}));
          connect(lF.y, BuckleyLeverett.u) annotation (points=[82.6,-30; 86,-30;
                86,68; -74,68; -74,56; -68.09,56],       style(color=74,
                rgbcolor={0,0,127}));
          connect(IC.y, BuckleyLeverett.u2) annotation (points=[-79.6,40;
                -68.09,40], style(color=74, rgbcolor={0,0,127}));
          connect(BCL.y, BuckleyLeverett.u3) annotation (points=[-79.6,30;
                -68.09,30], style(color=74, rgbcolor={0,0,127}));
          connect(BCR.y, BuckleyLeverett.u4) annotation (points=[-79.6,20; -74,
                20; -74,24; -68.09,24], style(color=74, rgbcolor={0,0,127}));
        end BuckleyLeverettEquation;
        annotation (Documentation(info="<html>
<p>
This package contains Buckley-Leverett equation solved with the Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end BuckleyLeverett;

      package EulerSystem
        package BaseModel
          model ShockWaveLDLR
            PDE.FiniteVolume.FVMIntegrator.FVIntegrator Density
              annotation (extent=[-60,60; -34,80]);
            PDE.FiniteVolume.FVMIntegrator.FVIntegrator Momentum
              annotation (extent=[-60,0; -34,20]);
            PDE.FiniteVolume.FVMIntegrator.FVIntegrator Energy
              annotation (extent=[-60,-60; -34,-40]);
            PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus
                                               annotation (extent=[-12,72; 10,80]);
            PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus
                                             annotation (extent=[-12,58; 10,66]);
            PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus
                                           annotation (extent=[26,72; 46,80]);
            PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus
                                         annotation (extent=[26,58; 46,66]);
            annotation (Diagram, Documentation(info="<html>
<p>
Base model for the Euler equations.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus1
              annotation (extent=[-12,12; 10,20]);
            PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus1
                                              annotation (extent=[-12,-2; 10,6]);
            PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus1
                                            annotation (extent=[26,12; 46,20]);
            PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus1
                                          annotation (extent=[26,-2; 46,6]);
            PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus2
              annotation (extent=[-12,-48; 10,-40]);
            PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus2
              annotation (extent=[-12,-62; 10,-54]);
            PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus2
                                            annotation (extent=[26,-48; 46,-40]);
            PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus2
                                          annotation (extent=[26,-62; 46,-54]);
            Modelica.Blocks.Math.Product product[worldModel1.n + 1]
              annotation (extent=[56,72; 64,80]);
            Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
              annotation (extent=[56,58; 64,66]);
            PDE.FiniteVolume.Fluxes.LaxFriedrichFlux.LF lF
                       annotation (extent=[74,62; 88,76]);
            Modelica.Blocks.Math.Division division[worldModel1.n + 1]
              annotation (extent=[26,40; 34,48]);
            Modelica.Blocks.Math.Division division1[worldModel1.n + 1]
              annotation (extent=[26,28; 34,36]);
            Modelica.Blocks.Math.Product product2[worldModel1.n + 1]
              annotation (extent=[56,12; 64,20]);
            Modelica.Blocks.Math.Product product3[worldModel1.n + 1]
              annotation (extent=[56,0; 64,8]);
            Modelica.Blocks.Math.Product product4[worldModel1.n + 1]
              annotation (extent=[-60,-80; -52,-72]);
            Modelica.Blocks.Math.Product product5[worldModel1.n + 1]
              annotation (extent=[-60,-96; -52,-88]);
            Modelica.Blocks.Sources.RealExpression const[worldModel1.n + 1](y=0.5)
              annotation (extent=[-74,-88; -68,-74]);
            Modelica.Blocks.Sources.RealExpression const1[worldModel1.n + 1](y=
                  0.5) annotation (extent=[-74,-104; -68,-90]);
            Modelica.Blocks.Math.Add add[worldModel1.n + 1](k2=-1)
              annotation (extent=[-40,-80; -32,-72]);
            Modelica.Blocks.Math.Add add1[worldModel1.n + 1](k2=-1)
              annotation (extent=[-40,-96; -32,-88]);
            Modelica.Blocks.Math.Product product6[worldModel1.n + 1]
              annotation (extent=[-10,-80; -2,-72]);
            Modelica.Blocks.Math.Product product7[worldModel1.n + 1]
              annotation (extent=[-10,-96; -2,-88]);
            Modelica.Blocks.Sources.RealExpression gammaMinusOne[worldModel1.n +
              1](y=0.4) annotation (extent=[-24,-90; -18,-78]);
            Modelica.Blocks.Math.Add add2[worldModel1.n + 1]
              annotation (extent=[74,12; 82,20]);
            Modelica.Blocks.Math.Add add3[worldModel1.n + 1]
              annotation (extent=[74,0; 82,8]);
            PDE.FiniteVolume.Fluxes.LaxFriedrichFlux.LF lF1
                        annotation (extent=[88,6; 98,16]);
            Modelica.Blocks.Math.Add add4[worldModel1.n + 1]
              annotation (extent=[60,-48; 68,-40]);
            Modelica.Blocks.Math.Add add5[worldModel1.n + 1]
              annotation (extent=[60,-62; 68,-54]);
            Modelica.Blocks.Math.Product product8[worldModel1.n + 1]
              annotation (extent=[76,-48; 84,-40]);
            Modelica.Blocks.Math.Product product9[worldModel1.n + 1]
              annotation (extent=[76,-62; 84,-54]);
            PDE.FiniteVolume.Fluxes.LaxFriedrichFlux.LF lF2
                        annotation (extent=[90,-54; 98,-46]);
            inner World.worldModel worldModel1(n=10)
              annotation (extent=[-20,90; 20,100]);
          equation
            connect(Density.y, lDLRminus_plus.u) annotation (points=[-32.7,76;
                  -14.2,76], style(color=74, rgbcolor={0,0,127}));
            connect(Density.y, lDLRplus_plus.u) annotation (points=[-32.7,76; -20,
                  76; -20,62; -14.2,62], style(color=74, rgbcolor={0,0,127}));
            connect(lDLRminus_plus.y, u_minus_plus.u1) annotation (points=[11.1,
                  76; 20,76; 20,73.6; 24,73.6], style(color=74, rgbcolor={0,0,127}));
            connect(lDLRplus_plus.y, u_plus_plus.u1) annotation (points=[11.1,62;
                  20,62; 20,59.6; 24,59.6], style(color=74, rgbcolor={0,0,127}));
            connect(Density.y, u_minus_plus.u) annotation (points=[-32.7,76; -20,
                  76; -20,82; 20,82; 20,78.4; 24,78.4], style(color=74, rgbcolor=
                    {0,0,127}));
            connect(Density.y, u_plus_plus.u) annotation (points=[-32.7,76; -20,
                  76; -20,68; 20,68; 20,64.4; 24,64.4], style(color=74, rgbcolor=
                    {0,0,127}));
            connect(Momentum.y, lDLRminus_plus1.u) annotation (points=[-32.7,16;
                  -14.2,16], style(color=74, rgbcolor={0,0,127}));
            connect(Momentum.y, lDLRplus_plus1.u) annotation (points=[-32.7,16;
                  -20,16; -20,2; -14.2,2], style(color=74, rgbcolor={0,0,127}));
            connect(Momentum.y, u_minus_plus1.u) annotation (points=[-32.7,16;
                  -20,16; -20,22; 20,22; 20,18.4; 24,18.4], style(color=74,
                  rgbcolor={0,0,127}));
            connect(lDLRminus_plus1.y, u_minus_plus1.u1) annotation (points=[11.1,
                  16; 20,16; 20,13.6; 24,13.6], style(color=74, rgbcolor={0,0,127}));
            connect(Momentum.y, u_plus_plus1.u) annotation (points=[-32.7,16; -20,
                  16; -20,8; 20,8; 20,4.4; 24,4.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(lDLRplus_plus1.y, u_plus_plus1.u1) annotation (points=[11.1,2;
                  20,2; 20,-0.4; 24,-0.4], style(color=74, rgbcolor={0,0,127}));
            connect(Energy.y, lDLRminus_plus2.u) annotation (points=[-32.7,-44;
                  -14.2,-44], style(color=74, rgbcolor={0,0,127}));
            connect(Energy.y, lDLRplus_plus2.u) annotation (points=[-32.7,-44;
                  -20,-44; -20,-58; -14.2,-58], style(color=74, rgbcolor={0,0,127}));
            connect(Energy.y, u_minus_plus2.u) annotation (points=[-32.7,-44; -20,
                  -44; -20,-38; 20,-38; 20,-41.6; 24,-41.6], style(color=74,
                  rgbcolor={0,0,127}));
            connect(Energy.y, u_plus_plus2.u) annotation (points=[-32.7,-44; -20,
                  -44; -20,-52; 20,-52; 20,-55.6; 24,-55.6], style(color=74,
                  rgbcolor={0,0,127}));
            connect(lDLRminus_plus2.y, u_minus_plus2.u1) annotation (points=[11.1,
                  -44; 20,-44; 20,-46.4; 24,-46.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(lDLRplus_plus2.y, u_plus_plus2.u1) annotation (points=[11.1,
                  -58; 20,-58; 20,-60.4; 24,-60.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(u_minus_plus.y, product.u1) annotation (points=[47,76; 50,76;
                  50,78.4; 55.2,78.4], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus.y, product1.u1) annotation (points=[47,62; 50,62;
                  50,64.4; 55.2,64.4], style(color=74, rgbcolor={0,0,127}));
            connect(product.y, lF.u) annotation (points=[64.4,76; 68,76; 68,74.6;
                  73.37,74.6], style(color=74, rgbcolor={0,0,127}));
            connect(product1.y, lF.u1) annotation (points=[64.4,62; 68,62; 68,
                  72.5; 73.37,72.5], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus.y, lF.u2) annotation (points=[47,62; 50,62; 50,70;
                  70,70; 70,65.5; 73.37,65.5], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus.y, lF.u3) annotation (points=[47,76; 50,76; 50,
                  56; 70,56; 70,63.4; 73.37,63.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(lF.y, Density.u) annotation (points=[88.7,69; 92,69; 92,86;
                  -68,86; -68,78; -61.43,78], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus1.y, division.u1) annotation (points=[47,16; 50,
                  16; 50,26; 20,26; 20,46.4; 25.2,46.4], style(color=74, rgbcolor=
                   {0,0,127}));
            connect(u_plus_plus1.y, division1.u1) annotation (points=[47,2; 50,2;
                  50,24; 22,24; 22,34.4; 25.2,34.4], style(color=74, rgbcolor={0,
                    0,127}));
            connect(u_minus_plus.y, division.u2) annotation (points=[47,76; 50,76;
                  50,52; 22,52; 22,41.6; 25.2,41.6], style(color=74, rgbcolor={0,
                    0,127}));
            connect(u_plus_plus.y, division1.u2) annotation (points=[47,62; 50,62;
                  50,52; 18,52; 18,29.6; 25.2,29.6], style(color=74, rgbcolor={0,
                    0,127}));
            connect(division.y, product.u2) annotation (points=[34.4,44; 52,44;
                  52,73.6; 55.2,73.6], style(color=74, rgbcolor={0,0,127}));
            connect(division1.y, product1.u2) annotation (points=[34.4,32; 55.2,
                  32; 55.2,59.6], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus1.y, product2.u2) annotation (points=[47,16; 52,
                  16; 52,13.6; 55.2,13.6], style(color=74, rgbcolor={0,0,127}));
            connect(division.y, product2.u1) annotation (points=[34.4,44; 52,44;
                  52,18.4; 55.2,18.4], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus1.y, product3.u2) annotation (points=[47,2; 51.1,2;
                  51.1,1.6; 55.2,1.6], style(color=74, rgbcolor={0,0,127}));
            connect(division1.y, product3.u1) annotation (points=[34.4,32; 52,32;
                  52,6.4; 55.2,6.4], style(color=74, rgbcolor={0,0,127}));
            connect(product2.y, product4.u1) annotation (points=[64.4,16; 70,16;
                  70,-8; -98,-8; -98,-73.6; -60.8,-73.6], style(color=74,
                  rgbcolor={0,0,127}));
            connect(const.y, product4.u2) annotation (points=[-67.7,-81; -63.85,
                  -81; -63.85,-78.4; -60.8,-78.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(const1.y, product5.u2) annotation (points=[-67.7,-97; -64.85,
                  -97; -64.85,-94.4; -60.8,-94.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(product3.y, product5.u1) annotation (points=[64.4,4; 68,4; 68,
                  -6; -96,-6; -96,-89.6; -60.8,-89.6], style(color=74, rgbcolor={
                    0,0,127}));
            connect(product4.y, add.u2) annotation (points=[-51.6,-76; -46,-76;
                  -46,-78.4; -40.8,-78.4], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus2.y, add.u1) annotation (points=[47,-44; 52,-44;
                  52,-68; -46,-68; -46,-73.6; -40.8,-73.6], style(color=74,
                  rgbcolor={0,0,127}));
            connect(product5.y, add1.u2) annotation (points=[-51.6,-92; -46,-92;
                  -46,-94.4; -40.8,-94.4], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus2.y, add1.u1) annotation (points=[47,-58; 50,-58;
                  50,-66; -48,-66; -48,-89.6; -40.8,-89.6], style(color=74,
                  rgbcolor={0,0,127}));
            connect(gammaMinusOne.y, product6.u2) annotation (points=[-17.7,-84;
                  -16,-84; -16,-78.4; -10.8,-78.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(gammaMinusOne.y, product7.u1) annotation (points=[-17.7,-84;
                  -16,-84; -16,-89.6; -10.8,-89.6], style(color=74, rgbcolor={0,0,
                    127}));
            connect(add.y, product6.u1) annotation (points=[-31.6,-76; -16,-76;
                  -16,-73.6; -10.8,-73.6], style(color=74, rgbcolor={0,0,127}));
            connect(add1.y, product7.u2) annotation (points=[-31.6,-92; -16,-92;
                  -16,-94.4; -10.8,-94.4], style(color=74, rgbcolor={0,0,127}));
            connect(product2.y, add2.u1) annotation (points=[64.4,16; 70,16; 70,
                  18.4; 73.2,18.4], style(color=74, rgbcolor={0,0,127}));
            connect(product3.y, add3.u1) annotation (points=[64.4,4; 68,4; 68,6.4;
                  73.2,6.4], style(color=74, rgbcolor={0,0,127}));
            connect(product6.y, add2.u2) annotation (points=[-1.6,-76; 56,-76;
                  56,-10; 68,-10; 68,13.6; 73.2,13.6],
                                                    style(color=74, rgbcolor={0,0,
                    127}));
            connect(product7.y, add3.u2) annotation (points=[-1.6,-92; 56,-92;
                  56,-10; 70,-10; 70,1.6; 73.2,1.6],
                                                  style(color=74, rgbcolor={0,0,
                    127}));
            connect(add2.y, lF1.u) annotation (points=[82.4,16; 84.975,16; 84.975,
                  15; 87.55,15], style(color=74, rgbcolor={0,0,127}));
            connect(add3.y, lF1.u1) annotation (points=[82.4,4; 84,4; 84,13.5;
                  87.55,13.5], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus1.y, lF1.u2) annotation (points=[47,2; 52,2; 52,-2;
                  86,-2; 86,8.5; 87.55,8.5], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus1.y, lF1.u3) annotation (points=[47,16; 50,16; 50,
                  -4; 87.55,-4; 87.55,7], style(color=74, rgbcolor={0,0,127}));
            connect(lF1.y, Momentum.u) annotation (points=[98.5,11; 100,11; 100,
                  24; -68,24; -68,18; -61.43,18], style(color=74, rgbcolor={0,0,
                    127}));
            connect(u_minus_plus2.y, add4.u1) annotation (points=[47,-44; 52,-44;
                  52,-41.6; 59.2,-41.6], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus2.y, add5.u1) annotation (points=[47,-58; 50,-58;
                  50,-55.6; 59.2,-55.6], style(color=74, rgbcolor={0,0,127}));
            connect(product6.y, add4.u2) annotation (points=[-1.6,-76; 54,-76;
                  54,-46.4; 59.2,-46.4],
                                      style(color=74, rgbcolor={0,0,127}));
            connect(product7.y, add5.u2) annotation (points=[-1.6,-92; 59.2,-92;
                  59.2,-60.4], style(color=74, rgbcolor={0,0,127}));
            connect(add4.y, product8.u2) annotation (points=[68.4,-44; 72,-44; 72,
                  -46.4; 75.2,-46.4], style(color=74, rgbcolor={0,0,127}));
            connect(add5.y, product9.u2) annotation (points=[68.4,-58; 72,-58; 72,
                  -60.4; 75.2,-60.4], style(color=74, rgbcolor={0,0,127}));
            connect(division.y, product8.u1) annotation (points=[34.4,44; 52,44;
                  52,-36; 72,-36; 72,-41.6; 75.2,-41.6], style(color=74, rgbcolor=
                   {0,0,127}));
            connect(division1.y, product9.u1) annotation (points=[34.4,32; 52,32;
                  52,-34; 75.2,-34; 75.2,-55.6], style(color=74, rgbcolor={0,0,
                    127}));
            connect(product8.y, lF2.u) annotation (points=[84.4,-44; 86,-44; 86,
                  -46.8; 89.64,-46.8], style(color=74, rgbcolor={0,0,127}));
            connect(product9.y, lF2.u1) annotation (points=[84.4,-58; 86,-58; 86,
                  -48; 89.64,-48], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus2.y, lF2.u2) annotation (points=[47,-58; 50,-58;
                  50,-52; 89.64,-52], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus2.y, lF2.u3) annotation (points=[47,-44; 52,-44;
                  52,-64; 88,-64; 88,-53.2; 89.64,-53.2], style(color=74,
                  rgbcolor={0,0,127}));
            connect(lF2.y, Energy.u) annotation (points=[98.4,-50; 100,-50; 100,
                  -32; -68,-32; -68,-42; -61.43,-42], style(color=74, rgbcolor={0,
                    0,127}));
          end ShockWaveLDLR;
          annotation (Documentation(info="<html>
<p>
This package contains the base model for the Euler equations.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end BaseModel;

        package SODProblem
          package LaxFriedrich
            model SOD
              extends
                PDE.FiniteVolume.Examples.EulerSystem.BaseModel.ShockWaveLDLR(
                                    worldModel1(n=10),
                Density(
                  vb=worldModel1.gcl + 1,
                  icb=worldModel1.gcl + 1,
                  ve=worldModel1.gcl + worldModel1.n - 1,
                  ice=worldModel1.gcl + worldModel1.n - 1,
                  bcl=0),
                Momentum(
                  vb=worldModel1.gcl + 2,
                  icb=worldModel1.gcl + 2,
                  bcr=0),
                Energy(
                  vb=worldModel1.gcl + 1,
                  icb=worldModel1.gcl + 1,
                  ve=worldModel1.gcl + worldModel1.n - 1,
                  ice=worldModel1.gcl + worldModel1.n - 1,
                  bcl=0),
                lF(alpha=0.2),
                lF1(alpha=0.2),
                lF2(alpha=0.2));
              Modelica.Blocks.Sources.RealExpression ICdensity[worldModel1.n](y=1.0)
                annotation (extent=[-80,64; -70,78]);
              Modelica.Blocks.Sources.RealExpression BCLdensity[worldModel1.gcl](y=0.0)
                         annotation (extent=[-80,52; -70,66]);
              Modelica.Blocks.Sources.RealExpression BCRdensity[worldModel1.gcr](y=
                    0.125) annotation (extent=[-80,40; -70,54]);
              annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Euler equations solved with Lax-Friedrichs flux</font></h3>
<p>
Implements the Euler system of equations
</p>
<p>
<img align=middle src=\"..\\Images\\euler.png\">
</p>

<p>
with Lax-Friedrichs flux. The initial conditions are
</p>

<img align=middle src=\"..\\Images\\sw4.png\">

<p>
The boundary conditions are
</p>

<img align=middle src=\"..\\Images\\sw6.png\">




</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
              Modelica.Blocks.Sources.RealExpression ICmomentum[worldModel1.n](y=
                    0.0) annotation (extent=[-90,16; -78,30]);
              Modelica.Blocks.Sources.RealExpression BCLmomentum[worldModel1.gcl](y=
                   0.0) annotation (extent=[-90,4; -78,18]);
              Modelica.Blocks.Sources.RealExpression BCRmomentum[worldModel1.gcr](y=
                   0.0) annotation (extent=[-90,-8; -78,6]);
              Modelica.Blocks.Sources.RealExpression ICenergy[worldModel1.n](y=2.5)
                annotation (extent=[-90,-56; -80,-42]);
              Modelica.Blocks.Sources.RealExpression BCLenergy[worldModel1.gcl](y=0.0)
                         annotation (extent=[-90,-66; -80,-52]);
              Modelica.Blocks.Sources.RealExpression BCRenergy[worldModel1.gcr](y=
                    0.25) annotation (extent=[-90,-78; -80,-62]);
            equation
              connect(ICdensity.y, Density.u2) annotation (points=[-69.5,71; -65.75,
                    71; -65.75,70; -61.43,70], style(color=74, rgbcolor={0,0,127}));
              connect(BCRdensity.y, Density.u4) annotation (points=[-69.5,47; -66,
                    47; -66,58; -64,58; -64,62; -61.43,62], style(color=74,
                    rgbcolor={0,0,127}));
              connect(ICmomentum.y, Momentum.u2) annotation (points=[-77.4,23; -70,
                    23; -70,16; -66,16; -66,10; -61.43,10], style(color=74,
                    rgbcolor={0,0,127}));
              connect(BCLmomentum.y, Momentum.u3) annotation (points=[-77.4,11; -70,
                    11; -70,5; -61.43,5], style(color=74, rgbcolor={0,0,127}));
              connect(ICenergy.y, Energy.u2) annotation (points=[-79.5,-49; -68,-49;
                    -68,-50; -61.43,-50], style(color=74, rgbcolor={0,0,127}));
              connect(BCRenergy.y, Energy.u4) annotation (points=[-79.5,-70; -66,
                    -70; -66,-58; -61.43,-58], style(color=74, rgbcolor={0,0,127}));
              connect(BCLenergy.y, Energy.u3) annotation (points=[-79.5,-59;
                    -70.75,-59; -70.75,-55; -61.43,-55], style(color=74, rgbcolor=
                     {0,0,127}));
              connect(BCRmomentum.y, Momentum.u4) annotation (points=[-77.4,-1;
                    -69.7,-1; -69.7,2; -61.43,2], style(color=74, rgbcolor={0,0,
                      127}));
              connect(BCLdensity.y, Density.u3) annotation (points=[-69.5,59;
                    -65.75,59; -65.75,65; -61.43,65], style(color=74, rgbcolor={0,
                      0,127}));
            end SOD;
            annotation (Documentation(info="<html>
<p>
This package contains Euler equations solved with the Finite Volume Methods by using the Lax-Friedrichs flux.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          end LaxFriedrich;

          package Roe

            model Euler
              LDLR.Reconstruction.Rec rec
                                     annotation (extent=[-40,-24; -26,-14]);
              LDLR.Reconstruction.Rec rec1
                                      annotation (extent=[-40,-54; -26,-44]);
              LDLR.Reconstruction.Rec rec2
                                      annotation (extent=[-40,-84; -26,-74]);
              annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Euler equations solved with Roe큦 flux</font></h3>
<p>
Implements the Euler system of equations
</p>
<p>
<img align=middle src=\"..\\Images\\euler.png\">
</p>

<p>
with Roe큦 flux. The initial conditions are
</p>

<img align=middle src=\"..\\Images\\sw4.png\">

<p>
The boundary conditions are
</p>

<img align=middle src=\"..\\Images\\sw6.png\">

<p>
This system has eigenvalues
</p>

<img align=middle src=\"..\\Images\\eulerroe1.png\">

<p>
and eigenvector matrix
</p>

<img align=middle src=\"..\\Images\\eulerroe2.png\">

<p>
that are computed at each time step at each interface by the corresponding blocks.
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
              PDE.FiniteVolume.Fluxes.Roe.Averages.Vaverage vaverage
                                         annotation (extent=[-72,80; -60,92]);
              PDE.FiniteVolume.Fluxes.Roe.Averages.Aaverage aaverage
                                         annotation (extent=[-72,64; -60,76]);
              PDE.FiniteVolume.Fluxes.Roe.Averages.Haverage haverage
                                         annotation (extent=[-72,48; -60,60]);
              PDE.FiniteVolume.Fluxes.Roe.Averages.Daverage daverage
                                         annotation (extent=[-72,32; -60,44]);
              inner World.worldModel worldModel1(m=3,
                n=10,
                deltat=0)
                annotation (extent=[-20,90; 20,100]);
              PDE.FiniteVolume.Fluxes.Roe.Lambda.Lambdas lambdas
                                     annotation (extent=[-42,70; -24,86]);
              PDE.FiniteVolume.Fluxes.Roe.Wave.Waves waves
                               annotation (extent=[-42,48; -24,64]);
              PDE.FiniteVolume.Fluxes.Roe.WaveStrength.a a
                               annotation (extent=[-42,18; -18,40]);
              PDE.FiniteVolume.Fluxes.Roe.DeltaU.Deltau deltau
                                   annotation (extent=[-72,20; -60,28]);
              PDE.FiniteVolume.Fluxes.Roe.DeltaU.Deltau deltau1
                                    annotation (extent=[-72,10; -60,18]);
              PDE.FiniteVolume.Fluxes.Roe.DeltaU.Deltau deltau2
                                    annotation (extent=[-72,0; -60,8]);
              PDE.FiniteVolume.Fluxes.Roe.FluxDifference.FluxDiff fluxDiff
                                               annotation (extent=[18,66; 40,86]);
              PDE.FiniteVolume.Fluxes.Roe.FluxDifference.FluxDiff fluxDiff1
                                                annotation (extent=[18,40; 40,60]);
              PDE.FiniteVolume.Fluxes.Roe.FluxDifference.FluxDiff fluxDiff2
                                                annotation (extent=[18,14; 40,34]);
              Modelica.Blocks.Math.Add3 add3_1[worldModel1.m,worldModel1.n + 1]
                annotation (extent=[60,52; 80,72]);
              Modelica.Blocks.Math.Add3 add3_2[worldModel1.m,worldModel1.n + 1]
                annotation (extent=[60,26; 80,46]);
              Modelica.Blocks.Math.Product product[worldModel1.n + 1]
                annotation (extent=[4,-40; 12,-32]);
              Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
                annotation (extent=[4,-60; 12,-52]);
              Modelica.Blocks.Math.Division vplus[worldModel1.n + 1]
                annotation (extent=[-14,-60; -6,-52]);
              Modelica.Blocks.Math.Division vminus[worldModel1.n + 1]
                annotation (extent=[-14,-40; -6,-32]);
              Modelica.Blocks.Math.Product product2[worldModel1.n + 1]
                annotation (extent=[24,-40; 32,-32]);
              Modelica.Blocks.Math.Product product3[worldModel1.n + 1]
                annotation (extent=[24,-60; 32,-52]);
              Modelica.Blocks.Math.Product product4[worldModel1.n + 1]
                annotation (extent=[42,-40; 50,-32]);
              Modelica.Blocks.Math.Product product5[worldModel1.n + 1]
                annotation (extent=[42,-60; 50,-52]);
              Modelica.Blocks.Sources.RealExpression const[worldModel1.n + 1](y=0.5)
                annotation (extent=[28,-50; 36,-40]);
              Modelica.Blocks.Math.Add add[worldModel1.n + 1](k2=-1)
                annotation (extent=[62,-40; 70,-32]);
              Modelica.Blocks.Math.Add add1[worldModel1.n + 1](k2=-1)
                annotation (extent=[62,-60; 70,-52]);
              Modelica.Blocks.Math.Product pminus[worldModel1.n + 1]
                annotation (extent=[80,-40; 88,-32]);
              Modelica.Blocks.Math.Product pplus[worldModel1.n + 1]
                annotation (extent=[80,-60; 88,-52]);
              Modelica.Blocks.Sources.RealExpression const1[worldModel1.n + 1](y=
                    0.4) annotation (extent=[62,-52; 70,-40]);
              Modelica.Blocks.Math.Add add2[worldModel1.n + 1]
                annotation (extent=[0,-74; 8,-66]);
              Modelica.Blocks.Math.Add add3[worldModel1.n + 1]
                annotation (extent=[0,-94; 8,-86]);
              Modelica.Blocks.Math.Division hminus[worldModel1.n + 1]
                annotation (extent=[62,-74; 70,-66]);
              Modelica.Blocks.Math.Division hplus[worldModel1.n + 1]
                annotation (extent=[62,-94; 70,-86]);
              PDE.FiniteVolume.Fluxes.Roe.IntegratorRoe.Integrator integrator(
                                                           bcl=0, bcr=0)
                annotation (extent=[-82,-94; -60,-74]);
              PDE.FiniteVolume.Fluxes.Roe.AverageFilter.Filter Density
                                           annotation (extent=[-84,-26; -60,-14]);
              PDE.FiniteVolume.Fluxes.Roe.AverageFilter.Filter Momentum(
                                                     row=2)
                annotation (extent=[-84,-48; -60,-36]);
              PDE.FiniteVolume.Fluxes.Roe.AverageFilter.Filter Energy(
                                                   row=3)
                annotation (extent=[-84,-70; -60,-58]);
              PDE.FiniteVolume.Examples.EulerSystem.SODProblem.Roe.IC iC
                    annotation (extent=[-100,-80; -94,-74]);
              PDE.FiniteVolume.Examples.EulerSystem.SODProblem.Roe.BCL bCL
                      annotation (extent=[-100,-90; -94,-84]);
              PDE.FiniteVolume.Examples.EulerSystem.SODProblem.Roe.BCR bCR
                      annotation (extent=[-100,-100; -94,-94]);
            equation
              connect(vaverage.y, lambdas.u) annotation (points=[-59.4,86; -52,86;
                    -52,82.88; -43.08,82.88], style(color=74, rgbcolor={0,0,127}));
              connect(aaverage.y, lambdas.u1) annotation (points=[-59.4,70; -48,70;
                    -48,73.12; -43.08,73.12], style(color=74, rgbcolor={0,0,127}));
              connect(vaverage.y, waves.u) annotation (points=[-59.4,86; -52,86;
                    -52,60.88; -43.17,60.88], style(color=74, rgbcolor={0,0,127}));
              connect(aaverage.y, waves.u1) annotation (points=[-59.4,70; -48,70;
                    -48,56.08; -43.17,56.08], style(color=74, rgbcolor={0,0,127}));
              connect(haverage.y, waves.u2) annotation (points=[-59.4,54; -56,54;
                    -56,51.2; -43.17,51.2], style(color=74, rgbcolor={0,0,127}));
              connect(aaverage.y, a.u1) annotation (points=[-59.4,70; -48,70; -48,
                    33.29; -43.08,33.29], style(color=74, rgbcolor={0,0,127}));
              connect(deltau.y, a.u2) annotation (points=[-59.4,24; -54,24; -54,
                    29.11; -43.08,29.11], style(color=74, rgbcolor={0,0,127}));
              connect(deltau1.y, a.u3) annotation (points=[-59.4,14; -50,14; -50,
                    24.6; -43.08,24.6], style(color=74, rgbcolor={0,0,127}));
              connect(deltau2.y, a.u4) annotation (points=[-59.4,4; -46,4; -46,
                    20.42; -43.08,20.42], style(color=74, rgbcolor={0,0,127}));
              connect(daverage.y, a.u) annotation (points=[-59.4,38; -51.24,38;
                    -51.24,37.8; -43.08,37.8], style(color=74, rgbcolor={0,0,127}));
              connect(lambdas.y, fluxDiff.u2) annotation (points=[-23.1,82.8;
                    -12,82.8; -12,70.2; 17.01,70.2],
                                                 style(color=74, rgbcolor={0,0,127}));
              connect(waves.y, fluxDiff.u1) annotation (points=[-23.1,60.8; -20,
                    60.8; -20,76.1; 17.01,76.1], style(color=74, rgbcolor={0,0,127}));
              connect(a.y, fluxDiff.u) annotation (points=[-16.8,35.6; -8,35.6;
                    -8,82.3; 17.01,82.3],
                                       style(color=74, rgbcolor={0,0,127}));
              connect(lambdas.y1, fluxDiff1.u2) annotation (points=[-23.1,78;
                    -14,78; -14,44.2; 17.01,44.2],
                                               style(color=74, rgbcolor={0,0,127}));
              connect(lambdas.y2, fluxDiff2.u2) annotation (points=[-23.1,73.2;
                    -4,73.2; -4,18.2; 17.01,18.2],
                                                style(color=74, rgbcolor={0,0,127}));
              connect(waves.y1y, fluxDiff1.u1) annotation (points=[-23.1,56; 0,
                    56; 0,50.1; 17.01,50.1],
                                         style(color=74, rgbcolor={0,0,127}));
              connect(waves.y2y, fluxDiff2.u1) annotation (points=[-23.1,51.2;
                    -2,51.2; -2,24.1; 17.01,24.1],
                                                style(color=74, rgbcolor={0,0,127}));
              connect(a.y1, fluxDiff1.u) annotation (points=[-16.8,29; 4,29; 4,
                    56.3; 17.01,56.3],
                                 style(color=74, rgbcolor={0,0,127}));
              connect(a.y2, fluxDiff2.u) annotation (points=[-16.8,22.4; 8,22.4;
                    8,30.3; 17.01,30.3],
                                       style(color=74, rgbcolor={0,0,127}));
              connect(fluxDiff.y, add3_1.u1) annotation (points=[41.1,82; 52,82; 52,
                    70; 58,70], style(color=74, rgbcolor={0,0,127}));
              connect(fluxDiff1.y, add3_1.u2) annotation (points=[41.1,56; 52,56;
                    52,62; 58,62], style(color=74, rgbcolor={0,0,127}));
              connect(fluxDiff2.y, add3_1.u3) annotation (points=[41.1,30; 52,30;
                    52,54; 58,54], style(color=74, rgbcolor={0,0,127}));
              connect(fluxDiff.y1, add3_2.u1) annotation (points=[41.1,70; 50,70;
                    50,44; 58,44], style(color=74, rgbcolor={0,0,127}));
              connect(fluxDiff1.y1, add3_2.u2) annotation (points=[41.1,44; 48,44;
                    48,36; 58,36], style(color=74, rgbcolor={0,0,127}));
              connect(fluxDiff2.y1, add3_2.u3) annotation (points=[41.1,18; 48,18;
                    48,28; 58,28], style(color=74, rgbcolor={0,0,127}));
              connect(rec1.y, vminus.u1) annotation (points=[-25.3,-46; -18,-46;
                    -18,-33.6; -14.8,-33.6], style(color=74, rgbcolor={0,0,127}));
              connect(rec1.y1, vplus.u1) annotation (points=[-25.3,-52; -18,-52;
                    -18,-53.6; -14.8,-53.6], style(color=74, rgbcolor={0,0,127}));
              connect(rec.y, vminus.u2) annotation (points=[-25.3,-16; -20,-16;
                    -20,-38.4; -14.8,-38.4],
                                         style(color=74, rgbcolor={0,0,127}));
              connect(rec.y1, vplus.u2) annotation (points=[-25.3,-22; -22,-22;
                    -22,-58.4; -14.8,-58.4],
                                         style(color=74, rgbcolor={0,0,127}));
              connect(vminus.y, product.u1) annotation (points=[-5.6,-36; -2,-36;
                    -2,-33.6; 3.2,-33.6], style(color=74, rgbcolor={0,0,127}));
              connect(vminus.y, product.u2) annotation (points=[-5.6,-36; -2,-36;
                    -2,-38.4; 3.2,-38.4], style(color=74, rgbcolor={0,0,127}));
              connect(vplus.y, product1.u1) annotation (points=[-5.6,-56; -2,-56;
                    -2,-53.6; 3.2,-53.6], style(color=74, rgbcolor={0,0,127}));
              connect(vplus.y, product1.u2) annotation (points=[-5.6,-56; -2,-56;
                    -2,-58.4; 3.2,-58.4], style(color=74, rgbcolor={0,0,127}));
              connect(product.y, product2.u2) annotation (points=[12.4,-36; 18,-36;
                    18,-38.4; 23.2,-38.4], style(color=74, rgbcolor={0,0,127}));
              connect(product1.y, product3.u2) annotation (points=[12.4,-56; 18,-56;
                    18,-58.4; 23.2,-58.4], style(color=74, rgbcolor={0,0,127}));
              connect(rec.y, product2.u1) annotation (points=[-25.3,-16; 18,-16;
                    18,-33.6; 23.2,-33.6],
                                        style(color=74, rgbcolor={0,0,127}));
              connect(rec.y1, product3.u1) annotation (points=[-25.3,-22; 16,
                    -22; 16,-53.6; 23.2,-53.6],
                                           style(color=74, rgbcolor={0,0,127}));
              connect(const.y, product4.u2) annotation (points=[36.4,-45; 38,-45;
                    38,-38.4; 41.2,-38.4], style(color=74, rgbcolor={0,0,127}));
              connect(const.y, product5.u1) annotation (points=[36.4,-45; 38,-45;
                    38,-53.6; 41.2,-53.6], style(color=74, rgbcolor={0,0,127}));
              connect(product2.y, product4.u1) annotation (points=[32.4,-36; 38,-36;
                    38,-33.6; 41.2,-33.6], style(color=74, rgbcolor={0,0,127}));
              connect(product3.y, product5.u2) annotation (points=[32.4,-56; 38,-56;
                    38,-58.4; 41.2,-58.4], style(color=74, rgbcolor={0,0,127}));
              connect(product4.y, add.u2) annotation (points=[50.4,-36; 58,-36; 58,
                    -38.4; 61.2,-38.4], style(color=74, rgbcolor={0,0,127}));
              connect(product5.y, add1.u2) annotation (points=[50.4,-56; 58,-56; 58,
                    -58.4; 61.2,-58.4], style(color=74, rgbcolor={0,0,127}));
              connect(rec2.y, add.u1) annotation (points=[-25.3,-76; 52,-76; 52,
                    -33.6; 61.2,-33.6], style(color=74, rgbcolor={0,0,127}));
              connect(rec2.y1, add1.u1) annotation (points=[-25.3,-82; 54,-82;
                    54,-53.6; 61.2,-53.6],
                                        style(color=74, rgbcolor={0,0,127}));
              connect(const1.y, pminus.u2) annotation (points=[70.4,-46; 76,-46; 76,
                    -38.4; 79.2,-38.4], style(color=74, rgbcolor={0,0,127}));
              connect(const1.y, pplus.u1) annotation (points=[70.4,-46; 76,-46; 76,
                    -53.6; 79.2,-53.6], style(color=74, rgbcolor={0,0,127}));
              connect(add.y, pminus.u1) annotation (points=[70.4,-36; 74,-36; 74,
                    -33.6; 79.2,-33.6], style(color=74, rgbcolor={0,0,127}));
              connect(add1.y, pplus.u2) annotation (points=[70.4,-56; 74,-56; 74,
                    -58.4; 79.2,-58.4], style(color=74, rgbcolor={0,0,127}));
              connect(rec2.y, add2.u1) annotation (points=[-25.3,-76; -20,-76;
                    -20,-67.6; -0.8,-67.6],
                                        style(color=74, rgbcolor={0,0,127}));
              connect(rec2.y1, add3.u1) annotation (points=[-25.3,-82; -20,-82;
                    -20,-87.6; -0.8,-87.6],
                                        style(color=74, rgbcolor={0,0,127}));
              connect(pminus.y, add2.u2) annotation (points=[88.4,-36; 92,-36;
                    92,-96; -8,-96; -8,-72.4; -0.8,-72.4],
                                                        style(color=74, rgbcolor={0,
                      0,127}));
              connect(pplus.y, add3.u2) annotation (points=[88.4,-56; 94,-56;
                    94,-98; -4,-98; -4,-92.4; -0.8,-92.4],
                                                        style(color=74, rgbcolor={0,
                      0,127}));
              connect(rec.y, vaverage.u) annotation (points=[-25.3,-16; -20,-16;
                    -20,-10; -98,-10; -98,90.8; -73.2,90.8], style(color=74,
                    rgbcolor={0,0,127}));
              connect(vminus.y, vaverage.u1) annotation (points=[-5.6,-36; -4,-36;
                    -4,-8; -96,-8; -96,88.4; -73.2,88.4], style(color=74, rgbcolor=
                      {0,0,127}));
              connect(rec.y1, vaverage.u2) annotation (points=[-25.3,-22; -22,
                    -22; -22,-6; -94,-6; -94,83.6; -73.2,83.6],
                                                           style(color=74, rgbcolor=
                     {0,0,127}));
              connect(vplus.y, vaverage.u3) annotation (points=[-5.6,-56; -2,-56;
                    -2,-4; -92,-4; -92,81.2; -73.2,81.2], style(color=74, rgbcolor=
                      {0,0,127}));
              connect(haverage.y, aaverage.u) annotation (points=[-59.4,54; -56,54;
                    -56,62; -80,62; -80,73.6; -73.2,73.6], style(color=74, rgbcolor=
                     {0,0,127}));
              connect(vaverage.y, aaverage.u1) annotation (points=[-59.4,86; -56,86;
                    -56,78; -82,78; -82,66.4; -73.2,66.4], style(color=74, rgbcolor=
                     {0,0,127}));
              connect(add2.y, hminus.u1) annotation (points=[8.4,-70; 14,-70; 14,
                    -78; 56,-78; 56,-67.6; 61.2,-67.6], style(color=74, rgbcolor={0,
                      0,127}));
              connect(add3.y, hplus.u1) annotation (points=[8.4,-90; 14,-90; 14,-84;
                    56,-84; 56,-87.6; 61.2,-87.6], style(color=74, rgbcolor={0,0,
                      127}));
              connect(rec.y, hminus.u2) annotation (points=[-25.3,-16; 18,-16;
                    18,4; 96,4; 96,-76; 58,-76; 58,-72.4; 61.2,-72.4],
                                                                 style(color=74,
                    rgbcolor={0,0,127}));
              connect(rec.y1, hplus.u2) annotation (points=[-25.3,-22; 16,-22;
                    16,6; 98,6; 98,-80; 58,-80; 58,-92.4; 61.2,-92.4],
                                                                 style(color=74,
                    rgbcolor={0,0,127}));
              connect(rec.y, haverage.u) annotation (points=[-25.3,-16; -20,-16;
                    -20,-10; -98,-10; -98,58.8; -73.2,58.8], style(color=74,
                    rgbcolor={0,0,127}));
              connect(hminus.y, haverage.u1) annotation (points=[70.4,-70; 100,-70;
                    100,8; 14,8; 14,-2; -90,-2; -90,56.4; -73.2,56.4], style(color=
                      74, rgbcolor={0,0,127}));
              connect(rec.y1, haverage.u2) annotation (points=[-25.3,-22; -22,
                    -22; -22,-6; -94,-6; -94,51.6; -73.2,51.6],
                                                           style(color=74, rgbcolor=
                     {0,0,127}));
              connect(hplus.y, haverage.u3) annotation (points=[70.4,-90; 100,-90;
                    100,10; 12,10; 12,-2; -88,-2; -88,49.2; -73.2,49.2], style(
                    color=74, rgbcolor={0,0,127}));
              connect(rec.y, daverage.u) annotation (points=[-25.3,-16; -22,-16;
                    -22,-6; -94,-6; -94,41.6; -73.2,41.6], style(color=74, rgbcolor=
                     {0,0,127}));
              connect(rec.y1, daverage.u1) annotation (points=[-25.3,-22; -20,
                    -22; -20,-10; -98,-10; -98,34.4; -73.2,34.4],
                                                             style(color=74,
                    rgbcolor={0,0,127}));
              connect(rec.y, deltau.u) annotation (points=[-25.3,-16; -22,-16;
                    -22,-6; -94,-6; -94,26.88; -72.72,26.88],
                                                          style(color=74, rgbcolor=
                      {0,0,127}));
              connect(rec.y1, deltau.u1) annotation (points=[-25.3,-22; -20,-22;
                    -20,-10; -98,-10; -98,21.12; -72.72,21.12], style(color=74,
                    rgbcolor={0,0,127}));
              connect(vminus.y, deltau1.u) annotation (points=[-5.6,-36; -4,-36; -4,
                    -2; -86,-2; -86,18; -72.72,18; -72.72,16.88], style(color=74,
                    rgbcolor={0,0,127}));
              connect(vplus.y, deltau1.u1) annotation (points=[-5.6,-56; -4,-56;
                    -4,-4; -84,-4; -84,11.12; -72.72,11.12],
                                                          style(color=74, rgbcolor=
                      {0,0,127}));
              connect(pminus.y, deltau2.u) annotation (points=[88.4,-36; 90,-36; 90,
                    2; -42,2; -42,-2; -82,-2; -82,6.88; -72.72,6.88], style(color=
                      74, rgbcolor={0,0,127}));
              connect(pplus.y, deltau2.u1) annotation (points=[88.4,-56; 92,-56; 92,
                    10; 10,10; 10,-2; -78,-2; -78,1.12; -72.72,1.12], style(color=
                      74, rgbcolor={0,0,127}));
              connect(Density.y, rec.u) annotation (points=[-58.8,-20; -50.1,-20;
                    -50.1,-19; -41.4,-19], style(color=74, rgbcolor={0,0,127}));
              connect(Momentum.y, rec1.u) annotation (points=[-58.8,-42; -50,-42;
                    -50,-49; -41.4,-49], style(color=74, rgbcolor={0,0,127}));
              connect(integrator.y, Density.u) annotation (points=[-58.9,-78; -54,
                    -78; -54,-72; -94,-72; -94,-20; -86.4,-20], style(color=74,
                    rgbcolor={0,0,127}));
              connect(integrator.y, Momentum.u) annotation (points=[-58.9,-78; -54,
                    -78; -54,-72; -94,-72; -94,-42; -86.4,-42], style(color=74,
                    rgbcolor={0,0,127}));
              connect(integrator.y, Energy.u) annotation (points=[-58.9,-78; -54,
                    -78; -54,-72; -94,-72; -94,-64; -86.4,-64], style(color=74,
                    rgbcolor={0,0,127}));
              connect(Energy.y, rec2.u) annotation (points=[-58.8,-64; -48,-64; -48,
                    -79; -41.4,-79], style(color=74, rgbcolor={0,0,127}));
              connect(add3_1.y, integrator.u) annotation (points=[81,62; 100,62;
                    100,-100; -88,-100; -88,-76.2; -82.55,-76.2], style(color=74,
                    rgbcolor={0,0,127}));
              connect(add3_2.y, integrator.u4) annotation (points=[81,36; 100,36;
                    100,-100; -86,-100; -86,-79; -82.55,-79], style(color=74,
                    rgbcolor={0,0,127}));
              connect(iC.y, integrator.u1) annotation (points=[-93.7,-77; -87.85,
                    -77; -87.85,-84.8; -82.55,-84.8], style(color=74, rgbcolor={0,0,
                      127}));
              connect(bCL.y, integrator.u2) annotation (points=[-93.7,-87; -87.85,
                    -87; -87.85,-89.2; -82.55,-89.2], style(color=74, rgbcolor={0,0,
                      127}));
              connect(bCR.y, integrator.u3) annotation (points=[-93.7,-97; -87.85,
                    -97; -87.85,-91.8; -82.55,-91.8], style(color=74, rgbcolor={0,0,
                      127}));
            end Euler;

            block IC
              extends Icons.BlockIcon;

            outer PDE.World.worldModel worldModel1;
            parameter Integer n = worldModel1.n;
            parameter Integer m = worldModel1.m;

            equation
            for j in 1:n loop
              y[1, j] = 1.0;
              y[2, j] = 0.0;
              y[3, j] = 2.5;
            end for;

            public
              Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n]
                annotation (extent=[100,-10; 120,10]);
              annotation (Icon(Text(
                    extent=[-44,34; 52,-28],
                    style(color=3, rgbcolor={0,0,255}),
                    string="ICroe")), Documentation(info="<html>
<p>
This block implements initial condition for the Euler equations.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            end IC;

            block BCL
              extends Icons.BlockIcon;

            outer PDE.World.worldModel worldModel1;
            parameter Integer gcl = worldModel1.gcl;

              Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.gcl]
                annotation (extent=[100,-10; 120,10]);
            equation
            for j in 1:gcl loop
              y[1, j] = 1.0;
              y[2, j] = 0.0;
              y[3, j] = 2.5;
            end for;

              annotation (Icon(Text(
                    extent=[-48,34; 62,-30],
                    style(color=3, rgbcolor={0,0,255}),
                    string="BCLroe")), Documentation(info="<html>
<p>
This block implements left boundary condition for the Euler equations.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            end BCL;

            block BCR
              extends Icons.BlockIcon;

            outer PDE.World.worldModel worldModel1;
            parameter Integer gcr = worldModel1.gcr;

            equation
            for j in 1:gcr loop
              y[1, j] = 0.125;
              y[2, j] = 0.0;
              y[3, j] = 0.25;
            end for;

            public
              Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.gcr]
                annotation (extent=[100,-10; 120,10]);
              annotation (Icon(Text(
                    extent=[-54,36; 64,-34],
                    style(color=3, rgbcolor={0,0,255}),
                    string="BCRroe")), Documentation(info="<html>
<p>
This block implements right boundary condition for the Euler equations.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            end BCR;
            annotation (Documentation(info="<html>
<p>
This package contains Euler equations solved with the Finite Volume Methods by using the Roe큦 flux.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          end Roe;
          annotation (Documentation(info="<html>
<p>
This package contains Euler equations solved with the Finite Volume Methods by using two numerical flux approximations: Lax-Friedrichs and Roe큦 fluxes.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end SODProblem;

        annotation (Documentation(info="<html>
<p>
This package contains Euler equations solved with the Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end EulerSystem;

      package Diffusion

        model DiffusionEquation
          PDE.FiniteVolume.FVMIntegrator.FVIntegrator Diffusion(
                                                  bcr=0,
            vb=3,
            icb=3,
            bcl=1,
            gcl=1,
            gcr=1)
            annotation (extent=[-66,0; -20,46]);
          annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Diffusion equation</font></h3>
<p>
Implements the diffusion equation
</p>

<img align=middle src=\"..\\Images\\d4.png\">


<p>
where sigma is a constant value. The initial condition is
</p>

<img align=middle src=\"..\\Images\\d1.png\">

<p>
and boundary conditions are
</p>

<img align=middle src=\"..\\Images\\d2.png\">

<p>
The analytical solution of this problem is implemented in DiffusionAnalytic block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Sources.RealExpression BCL[worldModel1.gcl](y=0.0)
                          annotation (extent=[-90,2; -82,16]);
          inner World.worldModel worldModel1(
            gcl=1,
            gcr=1,
            n=10)
            annotation (extent=[-20,90; 20,100]);
          PDE.FiniteVolume.Fluxes.DiffusionFlux.DiffusionFlux diffusionFlux(
                                      beta=0.01)
            annotation (extent=[6,28; 44,46]);
          Modelica.Blocks.Sources.RealExpression BCR[worldModel1.gcr]
            annotation (extent=[-90,-10; -82,4]);
          MOL.Examples.Diffusion.DIC1 dIC1_1
            annotation (extent=[-90,18; -80,28]);
        equation
          connect(Diffusion.y, diffusionFlux.u) annotation (points=[-17.7,36.8;
                -8,36.8; -8,37; 2.2,37],  style(color=74, rgbcolor={0,0,127}));
          connect(dIC1_1.y, Diffusion.u2) annotation (points=[-79.5,23; -73.75,23;
                -73.75,23; -68.53,23], style(color=74, rgbcolor={0,0,127}));
          connect(BCL.y, Diffusion.u3) annotation (points=[-81.6,9; -76,9; -76,
                11.5; -68.53,11.5], style(color=74, rgbcolor={0,0,127}));
          connect(BCR.y, Diffusion.u4) annotation (points=[-81.6,-3; -76,-3;
                -76,4.6; -68.53,4.6], style(color=74, rgbcolor={0,0,127}));
          connect(diffusionFlux.y, Diffusion.u) annotation (points=[45.9,37; 60,
                37; 60,60; -80,60; -80,41.4; -68.53,41.4], style(color=74,
                rgbcolor={0,0,127}));
        end DiffusionEquation;
        annotation (Documentation(info="<html>
<p>
This package contains diffusion equation solved with the Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Diffusion;

      package Advection

        package LaxFriedrich
          model AdvectionLF
            PDE.FiniteVolume.FVMIntegrator.FVIntegrator Advection(
                                                           bcr=0, bcl=0)
              annotation (extent=[-54,0; -24,28]);
            inner World.worldModel worldModel1(n=10)
              annotation (extent=[-20,88; 20,100]);
            annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Advection equation solved with Lax-Friedrichs flux</font></h3>
<p>
Implements the advection equation
</p>

<img align=middle src=\"..\\Images\\a1.png\">

<p>
where c is speed, with the Lax-Friedrichs numerical flux. The initial condition is
</p>

<img align=middle src=\"..\\Images\\a3.png\">

<p>
and boundary condition at the left is
</p>

<img align=middle src=\"..\\Images\\a4.png\">

<p>
The analytical solution of this problem is implemented in AdvectionAnalytic block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus
                                               annotation (extent=[-6,18; 16,28]);
            PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus
                                             annotation (extent=[-6,0; 16,10]);
            PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus
                                           annotation (extent=[34,18; 54,28]);
            PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus
                                         annotation (extent=[34,0; 54,10]);
            Modelica.Blocks.Math.Product product[worldModel1.n + 1]
              annotation (extent=[68,20; 76,28]);
            Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
              annotation (extent=[68,0; 76,8]);
            Modelica.Blocks.Sources.RealExpression const[worldModel1.n + 1](y=0.1)
              annotation (extent=[54,8; 62,20]);
            PDE.FiniteVolume.Fluxes.LaxFriedrichFlux.LF lF(
                                                alpha=0.2)
                       annotation (extent=[82,10; 92,20]);
            Modelica.Blocks.Sources.RealExpression BCR[worldModel1.gcr]
              annotation (extent=[-86,-18; -76,0]);
            MOL.Examples.Diffusion.DiffusionIC diffusionIC
              annotation (extent=[-86,20; -78,28]);
            Modelica.Blocks.Sources.RealExpression BCL[worldModel1.gcl](y=cos(-0.1
                  *time)) annotation (extent=[-90,-2; -72,16]);
          equation
            connect(Advection.y, lDLRminus_plus.u) annotation (points=[-22.5,
                  22.4; -16.25,22.4; -16.25,23; -8.2,23],
                                                    style(color=74, rgbcolor={0,0,
                    127}));
            connect(Advection.y, lDLRplus_plus.u) annotation (points=[-22.5,
                  22.4; -16,22.4; -16,5; -8.2,5],
                                            style(color=74, rgbcolor={0,0,127}));
            connect(lDLRminus_plus.y, u_minus_plus.u1) annotation (points=[17.1,23;
                  22,23; 22,20; 32,20],     style(color=74, rgbcolor={0,0,127}));
            connect(lDLRplus_plus.y, u_plus_plus.u1) annotation (points=[17.1,5;
                  22,5; 22,2; 32,2], style(color=74, rgbcolor={0,0,127}));
            connect(Advection.y, u_minus_plus.u) annotation (points=[-22.5,22.4;
                  -16,22.4; -16,32; 26,32; 26,26; 32,26], style(color=74,
                  rgbcolor={0,0,127}));
            connect(Advection.y, u_plus_plus.u) annotation (points=[-22.5,22.4;
                  -16,22.4; -16,14; 26,14; 26,8; 32,8], style(color=74, rgbcolor=
                    {0,0,127}));
            connect(const.y, product.u2) annotation (points=[62.4,14; 64,14; 64,
                  21.6; 67.2,21.6], style(color=74, rgbcolor={0,0,127}));
            connect(const.y, product1.u1) annotation (points=[62.4,14; 64,14; 64,
                  6.4; 67.2,6.4], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus.y, product1.u2) annotation (points=[55,5; 58,5;
                  58,1.6; 67.2,1.6], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus.y, product.u1) annotation (points=[55,23; 58,23;
                  58,26.4; 67.2,26.4], style(color=74, rgbcolor={0,0,127}));
            connect(product.y, lF.u) annotation (points=[76.4,24; 78,24; 78,19;
                  81.55,19], style(color=74, rgbcolor={0,0,127}));
            connect(product1.y, lF.u1) annotation (points=[76.4,4; 78,4; 78,17.5;
                  81.55,17.5], style(color=74, rgbcolor={0,0,127}));
            connect(u_plus_plus.y, lF.u2) annotation (points=[55,5; 58,5; 58,10;
                  80,10; 80,12.5; 81.55,12.5], style(color=74, rgbcolor={0,0,127}));
            connect(u_minus_plus.y, lF.u3) annotation (points=[55,23; 58,23; 58,
                  18; 66,18; 66,12; 76,12; 76,11; 81.55,11], style(color=74,
                  rgbcolor={0,0,127}));
            connect(lF.y, Advection.u) annotation (points=[92.5,15; 94,15; 94,
                  36; -62,36; -62,25.2; -55.65,25.2], style(color=74, rgbcolor={0,
                    0,127}));
            connect(diffusionIC.y, Advection.u2) annotation (points=[-77.6,24;
                  -66,24; -66,14; -55.65,14],
                                  style(color=74, rgbcolor={0,0,127}));
            connect(BCL.y, Advection.u3) annotation (points=[-71.1,7; -63.375,7;
                  -63.375,7; -55.65,7], style(color=74, rgbcolor={0,0,127}));
            connect(BCR.y, Advection.u4) annotation (points=[-75.5,-9; -66,-9;
                  -66,2.8; -55.65,2.8], style(color=74, rgbcolor={0,0,127}));
          end AdvectionLF;
          annotation (Documentation(info="<html>
<p>
This package contains advection equation solved with the Finite Volume Methods by using the Lax-Friedrichs flux.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end LaxFriedrich;

        package Upwind
          model AdvectionUpwind
            PDE.FiniteVolume.FVMIntegrator.FVIntegrator Advection(
                                                    bcl=0, bcr=0)
              annotation (extent=[-40,20; 0,60]);
            PDE.FiniteVolume.Fluxes.UpwindFlux.Upwind upwindNew
                                       annotation (extent=[30,36; 50,56]);
            annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Advection equation solved with Upwind flux</font></h3>
<p>
Implements the advection equation
</p>

<img align=middle src=\"..\\Images\\a1.png\">

<p>
where c is speed, with the Upwind numerical flux. The initial condition is
</p>

<img align=middle src=\"..\\Images\\a3.png\">

<p>
and boundary condition at the left is
</p>

<img align=middle src=\"..\\Images\\a4.png\">

<p>
The analytical solution of this problem is implemented in AdvectionAnalytic block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            Modelica.Blocks.Sources.RealExpression BCL[worldModel1.gcl](y=cos(-0.1
                  *time)) annotation (extent=[-80,20; -60,40]);
            Modelica.Blocks.Sources.RealExpression BCR[worldModel1.gcr]
              annotation (extent=[-76,6; -66,24]);
            Modelica.Blocks.Sources.RealExpression speed(y=0.1)
              annotation (extent=[6,32; 20,48]);
            MOL.Examples.Diffusion.DiffusionIC diffusionIC
              annotation (extent=[-76,42; -66,52]);
            inner World.worldModel worldModel1(n=10)
              annotation (extent=[-20,88; 20,100]);
          equation
            connect(Advection.y, upwindNew.u) annotation (points=[2,52; 28,52],
                style(color=74, rgbcolor={0,0,127}));
            connect(speed.y, upwindNew.u1)    annotation (points=[20.7,40; 28,
                  40],
                style(color=74, rgbcolor={0,0,127}));
            connect(BCL.y, Advection.u3) annotation (points=[-59,30; -42.2,30],
                style(color=74, rgbcolor={0,0,127}));
            connect(diffusionIC.y, Advection.u2) annotation (points=[-65.5,47;
                  -50.75,47; -50.75,40; -42.2,40], style(color=74, rgbcolor={0,0,
                    127}));
            connect(upwindNew.y, Advection.u) annotation (points=[51,46; 60,46;
                  60,64; -50,64; -50,56; -42.2,56], style(color=74, rgbcolor={0,0,
                    127}));
            connect(BCR.y, Advection.u4) annotation (points=[-65.5,15; -50,15;
                  -50,24; -42.2,24], style(color=74, rgbcolor={0,0,127}));
          end AdvectionUpwind;
          annotation (Documentation(info="<html>
<p>
This package contains advection equation solved with the Finite Volume Methods by using the upwind flux.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end Upwind;

        package LaxWendroff
          model AdvectionLW
            PDE.FiniteVolume.FVMIntegrator.FVIntegrator Advection(
                                                    bcl=0, bcr=0)
              annotation (extent=[-40,20; 0,60]);
            annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Advection equation solved with Lax-Wendroff flux</font></h3>
<p>
Implements the advection equation
</p>

<img align=middle src=\"..\\Images\\a1.png\">

<p>
where c is speed, with the Lax-Wendroff numerical flux. The initial condition is
</p>

<img align=middle src=\"..\\Images\\a3.png\">

<p>
and boundary condition at the left is
</p>

<img align=middle src=\"..\\Images\\a4.png\">

<p>
The analytical solution of this problem is implemented in AdvectionAnalytic block
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            Modelica.Blocks.Sources.RealExpression BCL[worldModel1.gcl](y=cos(-0.1
                  *time)) annotation (extent=[-80,20; -60,40]);
            Modelica.Blocks.Sources.RealExpression BCR[worldModel1.gcr]
              annotation (extent=[-80,-2; -60,18]);
            Modelica.Blocks.Sources.RealExpression velocity(y=0.1)
              annotation (extent=[6,22; 20,38]);
            MOL.Examples.Diffusion.DiffusionIC diffusionIC
              annotation (extent=[-76,42; -66,52]);
            inner World.worldModel worldModel1(deltat=0.2, n=40)
              annotation (extent=[-20,88; 20,100]);
            PDE.FiniteVolume.Fluxes.LaxWendroffFlux.LaxWendroff lW1_1
                             annotation (extent=[32,42; 44,54]);
          equation
            connect(BCR.y, Advection.u4) annotation (points=[-59,8; -50,8; -50,24;
                  -42.2,24], style(color=74, rgbcolor={0,0,127}));
            connect(BCL.y, Advection.u3) annotation (points=[-59,30; -42.2,30],
                style(color=74, rgbcolor={0,0,127}));
            connect(diffusionIC.y, Advection.u2) annotation (points=[-65.5,47;
                  -54.75,47; -54.75,40; -42.2,40], style(color=74, rgbcolor={0,0,
                    127}));
            connect(Advection.y, lW1_1.u) annotation (points=[2,52; 10.4,52; 10.4,
                  51.6; 30.8,51.6], style(color=74, rgbcolor={0,0,127}));
            connect(velocity.y, lW1_1.u1) annotation (points=[20.7,30; 26,30;
                  26,44.4; 30.8,44.4],
                                    style(color=74, rgbcolor={0,0,127}));
            connect(lW1_1.y, Advection.u) annotation (points=[44.6,48; 52,48; 52,
                  64; -50,64; -50,56; -42.2,56], style(color=74, rgbcolor={0,0,
                    127}));
          end AdvectionLW;
          annotation (Documentation(info="<html>
<p>
This package contains advection equation solved with the Finite Volume Methods by using the Lax-Wendroff flux.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end LaxWendroff;
        annotation (Documentation(info="<html>
<p>
This package contains advection equation solved with the Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Advection;

      annotation (Documentation(info="<html>
<p>
This package contains examples of partial differential equations solved with the Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end Examples;

    package FluxLimiter
      package FLIntegrator
        block FLIntegrator
          extends Icons.BlockIcon1;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Real beta = 1;
        parameter Integer vb = gcl+1 "|Unknowns| The left most unknown";
        parameter Integer ve = gcl + worldModel1.n
            "|Unknowns| The right most unknown";

        parameter Integer icb = gcl+1
            "|Initial Condition| Begin of the initial condition";
        parameter Integer ice = gcl + worldModel1.n
            "|Initial Condition| End of the initial condition";

        parameter Integer bcl = 1
            "|Boundary Conditions| Boundary condition at the left (0: no; 1: yes)";
        parameter Integer bcr = 1
            "|Boundary Conditions| Boundary condition at the right (0: no; 1: yes)";

        parameter Integer gcl = 2
            "|Boundary Conditions| Number of ghost cells at the left";
        parameter Integer gcr = 2
            "|Boundary Conditions| Number of ghost cells at the right";

        parameter Real delta_x = 1/n;
        parameter Real delta_t = 0.1;

        Real q[worldModel1.m, worldModel1.n + worldModel1.gcl + worldModel1.gcr];

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n + 1]
            annotation (extent=[-116,70; -100,90]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n + 1]
            annotation (extent=[-116,46; -100,66]);
          Modelica.Blocks.Interfaces.RealInput u2[worldModel1.m,worldModel1.n + 1]
            annotation (extent=[-116,16; -100,36]);
          Modelica.Blocks.Interfaces.RealInput u4[worldModel1.m,worldModel1.n]
            annotation (extent=[-116,-30; -100,-10]);
          Modelica.Blocks.Interfaces.RealInput u5[worldModel1.m,worldModel1.gcl]
            annotation (extent=[-116,-70; -100,-54]);
          Modelica.Blocks.Interfaces.RealInput u6[worldModel1.m,worldModel1.gcr]
            annotation (extent=[-116,-90; -100,-74]);
        equation
          y = q;

        for i in 1:m loop
          for j in 1:gcl loop
             q[i, j]= u5[i, j];
          end for;
        end for;

        for i in 1:m loop
          if bcl == 1 then
             q[i, gcl+1] = q[i, gcl];
          end if;
        end for;

        for i in 1:m loop
          for j in 1:gcr loop
             q[i, gcl+n+j] = u6[i, j];
          end for;
        end for;

        for i in 1:m loop
          if bcr == 1 then
             q[i, gcl+n] = q[i, gcl+n+1];
          end if;
        end for;

        for j in 1:m loop
          for i in vb:ve loop
            der(q[j, i]) = -(1/delta_x)*(u[j, i-gcl] + u1[j, i-gcl+1]) - (1/delta_x)*(u2[j, i-gcl+1] - u2[j, i-gcl]);
          end for;
        end for;

        initial equation

        for j in 1:m loop
          for i in icb:ice loop
            q[j, i] =  u4[j, i-gcl];
          end for;
        end for;

        public
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + worldModel1.gcl + worldModel1.gcr]
            annotation (extent=[100,50; 120,70]);
          annotation (Diagram, Icon(
              Text(
                extent=[-50,164; 48,94],
                style(color=3, rgbcolor={0,0,255}),
                string="%name"),
              Text(
                extent=[-108,88; -68,70],
                style(color=3, rgbcolor={0,0,255}),
                string="+"),
              Text(
                extent=[-102,68; -74,46],
                style(color=3, rgbcolor={0,0,255}),
                string="-"),
              Text(
                extent=[-92,34; -68,18],
                style(color=3, rgbcolor={0,0,255}),
                string="Flux"),
              Text(
                extent=[-100,-12; -68,-28],
                style(color=3, rgbcolor={0,0,255}),
                string="IC"),
              Text(
                extent=[-126,-54; -36,-70],
                style(color=3, rgbcolor={0,0,255}),
                string="gcl"),
              Text(
                extent=[-102,-74; -58,-90],
                style(color=3, rgbcolor={0,0,255}),
                string="gcr")),
            Documentation(info="<html>
<p>
The FLIntegrator takes as input the two fluctuation matrices, limited flux matrix, initial condition matrix and two boundary condition matrices. <br>
As output the block gives the average matrix <b>Q</b>. It implements the average update rule
</p>

<img align=middle src=\"..\\Images\\fl6.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end FLIntegrator;
        annotation (Documentation(info="<html>
<p>
This package contains integrator block that implements Finite Volume Methods by using flux limiters.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end FLIntegrator;

      package DeltaQ
        block deltaQ
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Integer p = 1;

        equation
        for i in 1:m loop
          for j in 1:n+gcl+gcr-1 loop
            y[i, j] = u[i, j+1] - u[i, j];
          end for;
        end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            worldModel1.gcl + worldModel1.gcr]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + worldModel1.gcl + worldModel1.gcr - 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(Bitmap(extent=[-46,46; 40,-44], name="Images/dq.PNG")),
                                  Diagram,
            Documentation(info="<html>
<p>
Takes the average matrix <b>Q</b> and computes the jumps
</p>

<img align=middle src=\"..\\Images\\fl9.png\">

<p>
giving matrix <b>deltaQ</b> as output. The entry delta_q[i, j] of the matrix contains the jump of the i-th unknown in the system at the interface j.
</p>
</pre>
<p><b>Release Notes: </b></p>

<ul>

</html>"));
        end deltaQ;
        annotation (Documentation(info="<html>
<p>
This package contains deltaQ block that computes the jumps at the cell boundaries.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end DeltaQ;

      package Teta
        block teta
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Integer p = 1;

          Modelica.Blocks.Interfaces.RealInput u1
            annotation (extent=[-140,-80; -100,-40]);
        equation
        // for i in 1:gcl loop
        //   y[i] = u[i];
        // end for;

        if u1 > 0 then
          for j in 1:n+1 loop
            if u[p, j+1] == 0 then
              y[p, j] = 0.0;
            else
              y[p, j] = u[p, j]/u[p, j+1];
            end if;
          end for;
        elseif u1 < 0 then
          for j in 1:n+1 loop
            if u[p, j+1] == 0 then
              y[p, j] = 0.0;
            else
              y[p, j] = u[p, j+2]/u[p, j+1];
            end if;
          end for;
        else
          for j in 1:n+1 loop
             y[p, j] = 0.0;
          end for;
        end if;

        // if u1 > 0 then
        //   for j in 1:n+1 loop
        //      y[p, j] = u[p, j]/u[p, j+1];
        //   end for;
        // elseif u1 < 0 then
        //   for j in 1:n+1 loop
        //      y[p, j] = u[p, j+2]/u[p, j+1];
        //   end for;
        // else
        //   for j in 1:n+1 loop
        //      y[p, j] = 0.0;
        //   end for;
        // end if;

        // for i in 1+gcl:gcl+n+1 loop
        //   if u1 > 0 then
        //     y[i-gcl] = (u[i-1] - u[i-2])/(u[i] - u[i-1]);
        //   else
        //     y[i-gcl] = (u[i+1] - u[i])/(u[i] - u[i-1]);
        //   end if;
        // end for;
        //
        // for i in gcl+n+1:gcl+n+gcr loop
        //   y[i] = u[i];
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            worldModel1.gcl + worldModel1.gcr - 1]
                                   annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Icon(
              Bitmap(extent=[-36,52; 32,-40], name="Images/teta.PNG"),
              Bitmap(extent=[-100,80; -62,40], name="Images/alpha.PNG"),
              Bitmap(extent=[-94,-38; -62,-88], name="Images/plambda.PNG")),
                                 Diagram,
            Documentation(info="<html>
<p>
Takes the <b>alpha</b> matrix and the j-th eigenvalue <b>lambdaj</b> as input and computes
</p>
<img align=middle src=\"..\\Images\\fl7.png\">
<p>
at each interface (i = 1,...,n+1). All other matrix entries (row # j) are set to zero. The resulting output is the matrix <b>teta</b>.
</p>
</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end teta;
        annotation (Documentation(info="<html>
<p>
This package contains tetha block that computes theta values at each cell interface.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Teta;

      package Wave

        block WaveP
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer m = worldModel1.m;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Integer p = 1;

        equation
        for j in 1:n+1 loop
          for i in 1:m loop
            y[i, j] = u[p, j]*u1[i, p];
          end for;
        end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            1] annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.m]
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1] annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(
              Text(
                extent=[-50,42; 64,-36],
                string="p-th Wave",
                style(color=0, rgbcolor={0,0,0})),
              Text(
                extent=[-102,-50; -58,-70],
                string="R",
                style(color=0, rgbcolor={0,0,0})),
              Bitmap(extent=[-98,82; -64,40], name="Images/alpha.PNG")),
            Documentation(info="<html>
<p>
Computes the matrix <b>j-th Waves</b> from either the matrix <b>alpha</b> or <b>Limitedalpha</b> and the eigenvalue matrix <b>R</b>.
</p>
<p>
<img align=middle src=\"..\\Images\\fl11.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\fl12.png\">
</p>

<p>
The <b>j-th Waves</b> matrix contains the j-th Wave of each interface.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end WaveP;
        annotation (Documentation(info="<html>
<p>
This package contains WaveP block that computes waves at each cell interface.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Wave;

      package Fluctuations
        block Fluctuation
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;

        equation
        if u > 0 then
          for j in 1:n+1 loop
            for i in 1:m loop
              y[i, j] = u*u1[i, j];
              y1[i, j] = 0;
            end for;
          end for;
        elseif u < 0 then
          for j in 1:n+1 loop
            for i in 1:m loop
              y1[i, j] = u*u1[i, j];
              y[i, j] = 0;
            end for;
          end for;
        else
        for j in 1:n+1 loop
            for i in 1:m loop
              y[i, j] = 0;
              y1[i, j] = 0;
            end for;
        end for;
        end if;

        public
          Modelica.Blocks.Interfaces.RealInput u
            annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[-140,-80; -100,-40]);
          annotation (Diagram, Icon(
              Text(
                extent=[48,56; 92,26],
                string="+",
                style(color=0, rgbcolor={0,0,0})),
              Text(
                extent=[54,-22; 88,-54],
                string="-",
                style(color=0, rgbcolor={0,0,0})),
              Text(
                extent=[-86,-34; -20,-82],
                string="p-th Waves",
                style(color=0, rgbcolor={0,0,0})),
              Text(
                extent=[-44,20; 42,-18],
                string="Fluctuation",
                style(color=0, rgbcolor={0,0,0})),
              Bitmap(extent=[-96,82; -54,36], name="Images/lambda.png")),
            Documentation(info="<html>
<p>
Computes two fluctuation matrices <b>A+</b> and <b>A-</b> from <b>j-th Waves</b> matrix and j-th eigenvalue <b>lambda</b>.
</p>
<p>
<img align=middle src=\"..\\Images\\fl13.png\">
</p>
<p>
<img align=middle src=\"..\\Images\\fl14.png\">
</p>
</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,30; 120,50]);
          Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,-50; 120,-30]);
        end Fluctuation;
        annotation (Documentation(info="<html>
<p>
This package contains Fluctuation block that computes positive and negative fluctuations at each cell interface.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Fluctuations;

      package Alpha
        block Alpha
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Integer p = 1;

        equation
          for j in 1:n+1 loop
            for i in 1:m loop
              y[i, j] = u[i, j+gcl-1];
            end for;
          end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            worldModel1.gcl + worldModel1.gcr - 1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1] annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Documentation(info="<html>
<p>
This block just filter the input matrix <b>alpha</b>, with boundary conditions, to the output matrix <b>alpha</b> without boundary conditions.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"),  Icon(Bitmap(extent=[-46,44; 44,-38], name="Images/alpha.PNG")));
        end Alpha;
        annotation (Documentation(info="<html>
<p>
This package contains Alpha block that filters the alpha matrix.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Alpha;

      package Limitedalpha
        block LimitedAlpha
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Integer p = 1;

        equation
        for j in 1:n+1 loop
          for i in 1:m loop
            y[i, j] = u[i, j]*u1[i, j];
            //y[i, j] = u[i, j]*u1[i, j+gcl-1];
          end for;
        end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            1]
            annotation (extent=[-140,40; -100,80]);
          annotation (Diagram, Icon(
              Bitmap(extent=[-56,52; 52,-36], name="Images/atilde.png"),
              Bitmap(extent=[-94,-36; -54,-84], name="Images/alpha.PNG"),
              Bitmap(extent=[-90,86; -46,34], name="Images/fteta.PNG")),
            Documentation(info="<html>
<p>
Computes the <b>limited alpha</b> matrix from matrices <b>f(teta)</b> and <b>alpha</b>.
</p>
<img align=middle src=\"..\\Images\\fl8.png\">
<p>
All other matrix entries (row # j) are set to zero.
</p>
</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[-140,-80; -100,-40]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,-10; 120,10]);
        end LimitedAlpha;
        annotation (Documentation(info="<html>
<p>
This package contains LimitedAlpha block that computes the limithed alpha values at each cell interface.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Limitedalpha;

      package Riemann
        block Riemann
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer m = worldModel1.m;
        parameter Integer n = worldModel1.n;

        Real alpha[worldModel1.m, worldModel1.n + worldModel1.gcl + worldModel1.gcr - 1];
        equation

          u = u1*alpha;
          y = alpha;
        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            worldModel1.gcl + worldModel1.gcr - 1]
               annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.m]
            annotation (extent=[-140,-80; -100,-40]);
          annotation (Diagram, Icon(
              Text(
                extent=[-108,-50; -56,-70],
                string="R",
                style(color=0, rgbcolor={0,0,0})),
              Text(
                extent=[-46,48; 48,-44],
                style(color=0, rgbcolor={0,0,0}),
                string="Riemann"),
              Bitmap(extent=[-94,82; -50,38], name="Images/dq.PNG")),
            Documentation(info="<html>
<p>
Takes the <b>deltaQ</b> matrix and the eigenvector matrix <b>R</b> and solves the Riemann problem deltaQ = R*alpha <br>
by finding the <b>alpha</b> matrix and giving it as output.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + worldModel1.gcl + worldModel1.gcr - 1]
                  annotation (extent=[100,-10; 120,10]);
        end Riemann;

        annotation (Documentation(info="<html>
<p>
This package contains Riemann block that solves the Riemann problem.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Riemann;

      package FluxLimited
        block FluxLimited
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        parameter Real deltat = worldModel1.deltat;
        parameter Real deltax = 1/n;
        parameter Integer p = 1;

        equation
        for j in 1:n+1 loop
          for i in 1:m loop
            y[i, j] = 0.5*abs(u)*(1-(deltat/deltax)*abs(u))*u1[i, j];
          end for;
        end for;

        public
          Modelica.Blocks.Interfaces.RealInput u
            annotation (extent=[-140,40; -100,80]);
          Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n
             + 1] annotation (extent=[-140,-80; -100,-40]);
          annotation (Diagram, Icon(
              Text(
                extent=[-40,34; 44,-32],
                string="Flux Limited",
                style(color=0, rgbcolor={0,0,0})),
              Text(
                extent=[-86,-44; -28,-76],
                string="p-th Wave",
                style(color=0, rgbcolor={0,0,0})),
              Bitmap(extent=[-98,74; -50,42], name="Images/plambda.PNG")),
            Documentation(info="<html>
<p>
Computes the j-th flux matrix <b>FluxLimited</b> at each interface from the j-th eigenvalue <b>j-th lambda</b> and the j-th wave matrix <b>j-th Waves</b>.
</p>
<img align=middle src=\"..\\Images\\fl10.png\">
<p>
Only one term in the above sum is computed. To obtain the complete sum, <b>m</b> FluxLimited blocks are needed. At the end all FluxLimited blocks must be summed.
</p>
</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1] annotation (extent=[100,-10; 120,10]);
        end FluxLimited;
        annotation (Documentation(info="<html>
<p>
This package contains FluxLimited block that computes fluxes at each cell interface.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end FluxLimited;

      package LinearMethods
        block Upwind
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;

        equation
          y = 0.0*u;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,-20; 140,20]);
          annotation (Diagram, Icon(Text(
                extent=[-50,30; 48,-26],
                style(color=3, rgbcolor={0,0,255}),
                string="Upwind")),
            Documentation(info="<html>
<p>
The upwind method takes as input the <b>teta</b> matrix and gives as output matrix
</p>

<img align=middle src=\"..\\Images\\fl1.png\">

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end Upwind;

        block LaxWendroff
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        outer parameter Integer p;

        equation
        for j in 1:n+1 loop
          for i in 1:p-1 loop
            y[i, j] = 0.0;
          end for;
        end for;

        for j in 1:n+1 loop
          y[p, j] = 1.0;
        end for;

        for j in 1:n+1 loop
          for i in p+1:m loop
            y[i, j] = 0.0;
          end for;
        end for;

        // for j in 1:n+1 loop
        //   for i in 1:m loop
        //     y[i, j] = 1.0;
        //   end for;
        // end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,-20; 140,20]);
          annotation (Diagram, Icon(Text(
                extent=[-68,54; 70,-50],
                style(color=3, rgbcolor={0,0,255}),
                string="LaxWendroff")),
            Documentation(info="<html>
<p>
The Lax Wendroff Method takes as input the <b>teta</b> matrix and gives as output the matrix
</p>

<img align=middle src=\"..\\Images\\fl2.png\">

<p>
which has ones in the j-th row and zeros elsewhere.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end LaxWendroff;

        block BeamWarming
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        outer parameter Integer p;

        equation
        for j in 1:n+1 loop
          for i in 1:p-1 loop
            y[i, j] = 0.0;
          end for;
        end for;

        for j in 1:n+1 loop
          y[p, j] = u[p, j];
        end for;

        for j in 1:n+1 loop
          for i in p+1:m loop
            y[i, j] = 0.0;
          end for;
        end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,-20; 140,20]);
          annotation (Diagram, Icon(Text(
                extent=[-82,64; 82,-60],
                style(color=3, rgbcolor={0,0,255}),
                string="BeamWarming")),
            Documentation(info="<html>
<p>
The Beam Warming Method takes as input the <b>teta</b> matrix and gives as output matrix
</p>

<img align=middle src=\"..\\Images\\fl3.png\">

<p>
which has in the j-th row the entries of the p-th row of the theta matrix and zeros in all other entries.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end BeamWarming;

        block Fromm
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        outer parameter Integer p;

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1]
            annotation (extent=[100,-20; 140,20]);
        equation

        for j in 1:n+1 loop
          for i in 1:p-1 loop
            y[i, j] = 0.0;
          end for;
        end for;

        for j in 1:n+1 loop
          y[p, j] = 0.5*(1+u[p, j]);
        end for;

        for j in 1:n+1 loop
          for i in p+1:m loop
            y[i, j] = 0.0;
          end for;
        end for;

          annotation (Diagram, Icon(Text(
                extent=[-48,28; 48,-24],
                style(color=3, rgbcolor={0,0,255}),
                string="Fromm")),
            Documentation(info="<html>
<p>
The Fromm Method takes as input the <b>teta</b> matrix and gives as output matrix <b>f(teta)</b> which has in the j-th row
</p>

<img align=middle src=\"..\\Images\\fl4.png\">

<p>
and zeros in all other rows.
</p>
</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end Fromm;
        annotation (Documentation(info="<html>
<p>
This package offers some linear methods for solving linear constant-coefficient hyperbolic system of equations.
</p>

</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end LinearMethods;

      package HighResolutionMethods

        block vanLeer
          extends Icons.BlockIcon;

        outer PDE.World.worldModel worldModel1;
        parameter Integer n = worldModel1.n;
        parameter Integer m = worldModel1.m;
        parameter Integer gcl = worldModel1.gcl;
        parameter Integer gcr = worldModel1.gcr;
        outer parameter Integer p;

        equation
        for j in 1:n+1 loop
          for i in 1:p-1 loop
            y[i, j] = 0.0;
          end for;
        end for;

        for j in 1:n+1 loop
          y[p, j] = (u[p, j] + abs(u[p, j]))/(1 + abs(u[p, j]));
        end for;

        for j in 1:n+1 loop
          for i in p+1:m loop
            y[i, j] = 0.0;
          end for;
        end for;

        public
          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
            1] annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
             + 1] annotation (extent=[100,-20; 140,20]);
          annotation (Diagram, Icon(Text(
                extent=[-54,38; 54,-34],
                style(color=3, rgbcolor={0,0,255}),
                string="van Leer")),
            Documentation(info="<html>
<p>
The van Leer Method takes as input the <b>teta</b> matrix and gives as output the <b>f(teta)</b> matrix which has in the j-th row
</p>

<img align=middle src=\"..\\Images\\fl5.png\">

<p>
and zeros in all other rows.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end vanLeer;
        annotation (Documentation(info="<html>
<p>
This package offers a high resolution method, van Leer method, for solving linear constant-coefficient hyperbolic system of equations.
</p>

</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end HighResolutionMethods;

      package FluxLimiterSolver
        block FluxLimiterSolver
          extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer method = worldModel1.fls;
          inner parameter Integer p = 1;

          Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n + 1]
            annotation (extent=[-140,-20; -100,20]);
          Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n + 1]
            annotation (extent=[100,-10; 120,10]);
          annotation (Diagram, Icon(Text(
                extent=[-44,34; 44,-30],
                string="FLSolver",
                style(color=0, rgbcolor={0,0,0}))),
            Documentation(info="<html>
<p>
The Flux Limiter Solver contains all the linear and high resolution methods which can be selected through the parameter variable <b>method</b> <br>
</p>

<p>
method = 1:    <b>Upwind</b> <br>
method = 2:    <b>Lax-Wendroff</b> <br>
method = 3:    <b>Beam Warming</b> <br>
method = 4:    <b>Fromm</b> <br>
method = 5:    <b>van Leer</b> <br>
</p>
</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          LinearMethods.Upwind upwind if method == 1
            annotation (extent=[-20,52; 20,74]);
          LinearMethods.LaxWendroff laxWendroff if method == 2
            annotation (extent=[-20,20; 20,44]);
          LinearMethods.BeamWarming beamWarming if method == 3
            annotation (extent=[-20,-12; 20,12]);
          LinearMethods.Fromm fromm if method == 4
            annotation (extent=[-20,-44; 20,-20]);
          PDE.FiniteVolume.FluxLimiter.HighResolutionMethods.vanLeer vanLeerNew
            if                                           method == 5
            annotation (extent=[-20,-76; 20,-52]);
        equation
          connect(u, upwind.u) annotation (points=[-120,0; -40,0; -40,63; -24,
                63], style(color=74, rgbcolor={0,0,127}));
          connect(u, laxWendroff.u) annotation (points=[-120,0; -40,0; -40,32;
                -24,32],   style(color=74, rgbcolor={0,0,127}));
          connect(u, beamWarming.u) annotation (points=[-120,0; -24,0], style(
                color=74, rgbcolor={0,0,127}));
          connect(u, fromm.u) annotation (points=[-120,0; -40,0; -40,-32; -24,
                -32], style(color=74, rgbcolor={0,0,127}));
          connect(u, vanLeerNew.u) annotation (points=[-120,0; -40,0; -40,-64;
                -24,-64], style(color=74, rgbcolor={0,0,127}));
          connect(upwind.y, y) annotation (points=[24,63; 40,63; 40,0; 110,0],
              style(color=74, rgbcolor={0,0,127}));
          connect(laxWendroff.y, y) annotation (points=[24,32; 40,32; 40,0; 110,
                0], style(color=74, rgbcolor={0,0,127}));
          connect(beamWarming.y, y) annotation (points=[24,0; 110,0], style(
                color=74, rgbcolor={0,0,127}));
          connect(fromm.y, y) annotation (points=[24,-32; 40,-32; 40,0; 110,0],
              style(color=74, rgbcolor={0,0,127}));
          connect(vanLeerNew.y, y) annotation (points=[24,-64; 40,-64; 40,0;
                110,0], style(color=74, rgbcolor={0,0,127}));
        end FluxLimiterSolver;
        annotation (Documentation(info="<html>
<p>
This package contains FluxLimiterSolver block that offers linear and high resolution methods for solving linear constant-coefficient hyperbolic system of equations.
</p>

</pre>
<p><b>Release Notes: </b></p>
<ul>
</html>"));
      end FluxLimiterSolver;

      package Examples
        package Acoustic
          model Acoustics
            annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Acoustics equations</font></h3>
<p>
The Acoustics example is a system of three equations
</p>

<img align=middle src=\"..\\Images\\ac1.png\">

<p>
where p is the gas pressure, v its velocity, phi the concentration of the contaminant and K<sub>0</sub>, rho<sub>0</sub>, v<sub>0</sub> are constant values.
The eigenvalues
</p>

<img align=middle src=\"..\\Images\\ac2.png\">
<img align=middle src=\"..\\Images\\ac3.png\">
<img align=middle src=\"..\\Images\\ac4.png\">

<p>
where c<sub>0</sub> is a constant value, and the eigenvector matrix
</p>

<img align=middle src=\"..\\Images\\ac5.png\">

<p>
are given explicitly in this example.
In general, the user must compute the eigenvector matrix and eigenvalues by using other Modelica packages.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            Riemann.Riemann riemann annotation (extent=[-70,66; -58,78]);
            Alpha.Alpha alpha annotation (extent=[-40,66; -28,78]);
            Wave.WaveP waveP annotation (extent=[-16,82; -4,94]);
            Wave.WaveP waveP1(p=2) annotation (extent=[-16,66; -4,78]);
            Wave.WaveP waveP2(p=3) annotation (extent=[-16,50; -4,62]);
            Fluctuations.Fluctuation fluctuation
              annotation (extent=[20,82; 32,94]);
            Fluctuations.Fluctuation fluctuation1
              annotation (extent=[20,66; 32,78]);
            Fluctuations.Fluctuation fluctuation2
              annotation (extent=[20,50; 32,62]);
            Modelica.Blocks.Math.Add3 add3_1[worldModel1.m,worldModel1.n + 1]
              annotation (extent=[50,74; 62,86]);
            Modelica.Blocks.Math.Add3 add3_2[worldModel1.m,worldModel1.n + 1]
              annotation (extent=[50,54; 62,66]);
            Teta.teta teta annotation (extent=[-60,0; -50,10]);
            Teta.teta teta1(p=2) annotation (extent=[-60,-16; -50,-6]);
            Teta.teta teta2(p=3) annotation (extent=[-60,-32; -50,-22]);
            Limitedalpha.LimitedAlpha limitedAlpha
              annotation (extent=[-20,0; -10,10]);
            Limitedalpha.LimitedAlpha limitedAlpha1(p=2)
              annotation (extent=[-20,-16; -10,-6]);
            Limitedalpha.LimitedAlpha limitedAlpha2(p=3)
              annotation (extent=[-20,-32; -10,-22]);
            Wave.WaveP waveP3 annotation (extent=[0,0; 10,10]);
            Wave.WaveP waveP4(p=2) annotation (extent=[0,-16; 10,-6]);
            Wave.WaveP waveP5(p=3) annotation (extent=[0,-32; 10,-22]);
            FluxLimited.FluxLimited fluxLimited annotation (extent=[20,0; 30,10]);
            FluxLimited.FluxLimited fluxLimited1(p=2)
              annotation (extent=[20,-16; 30,-6]);
            FluxLimited.FluxLimited fluxLimited2(p=3)
              annotation (extent=[20,-32; 30,-22]);
            Modelica.Blocks.Math.Add3 add3_3[worldModel1.m,worldModel1.n + 1]
              annotation (extent=[42,-18; 56,-4]);
            Modelica.Blocks.Sources.RealExpression lambda1(y=1.5)
              annotation (extent=[-96,34; -88,52]);
            Modelica.Blocks.Sources.RealExpression lambda2(y=2.0)
              annotation (extent=[-96,22; -88,40]);
            Modelica.Blocks.Sources.RealExpression lambda3(y=2.5)
              annotation (extent=[-96,10; -88,28]);
            Modelica.Blocks.Sources.RealExpression R[worldModel1.m,worldModel1.m](
               y={{-0.5,0,0.5},{1,0,1},{0,1,0}})
              annotation (extent=[-96,48; -80,66]);
            DeltaQ.deltaQ deltaQ annotation (extent=[-88,70; -80,78]);
            FLIntegrator.FLIntegrator Acoustics(bcl=0, bcr=0)
              annotation (extent=[34,-88; 86,-48]);
            inner World.worldModel worldModel1(
              m=3,
              deltat=0.1,
              n=10,
              fls=1)
              annotation (extent=[-100,-100; -60,-80]);
            FluxLimiterSolver.FluxLimiterSolver fluxLimiterSolver
              annotation (extent=[-40,0; -30,10]);
            FluxLimiterSolver.FluxLimiterSolver fluxLimiterSolver1(p=2)
              annotation (extent=[-40,-16; -30,-6]);
            FluxLimiterSolver.FluxLimiterSolver fluxLimiterSolver2(p=3)
              annotation (extent=[-40,-32; -30,-22]);
            Modelica.Blocks.Sources.RealExpression IC[worldModel1.m,worldModel1.
              n](y=2.0) annotation (extent=[12,-80; 20,-64]);
            Modelica.Blocks.Sources.RealExpression BCL[worldModel1.m,
              worldModel1.gcl](y=1.0) annotation (extent=[12,-92; 20,-76]);
            Modelica.Blocks.Sources.RealExpression BCR[worldModel1.m,
              worldModel1.gcr](y=1.0) annotation (extent=[12,-104; 20,-88]);
          equation
            connect(alpha.y, waveP.u) annotation (points=[-27.4,72; -24,72; -24,
                  91.6; -17.2,91.6], style(color=74, rgbcolor={0,0,127}));
            connect(alpha.y, waveP1.u) annotation (points=[-27.4,72; -24,72; -24,
                  75.6; -17.2,75.6], style(color=74, rgbcolor={0,0,127}));
            connect(alpha.y, waveP2.u) annotation (points=[-27.4,72; -24,72; -24,
                  59.6; -17.2,59.6], style(color=74, rgbcolor={0,0,127}));
            connect(waveP.y, fluctuation.u1) annotation (points=[-3.4,88; 8,88;
                  8,84.4; 18.8,84.4],
                                    style(color=74, rgbcolor={0,0,127}));
            connect(waveP1.y, fluctuation1.u1) annotation (points=[-3.4,72; 8,
                  72; 8,68.4; 18.8,68.4],
                                      style(color=74, rgbcolor={0,0,127}));
            connect(waveP2.y, fluctuation2.u1) annotation (points=[-3.4,56; 8,
                  56; 8,52.4; 18.8,52.4],
                                      style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation.y, add3_1.u1) annotation (points=[32.6,90.4; 40,
                  90.4; 40,84.8; 48.8,84.8], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1.y, add3_1.u2) annotation (points=[32.6,74.4; 40,
                  74.4; 40,80; 48.8,80], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation2.y, add3_1.u3) annotation (points=[32.6,58.4; 40,
                  58.4; 40,75.2; 48.8,75.2], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation.y1, add3_2.u1) annotation (points=[32.6,85.6; 38,
                  85.6; 38,64.8; 48.8,64.8], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1.y1, add3_2.u2) annotation (points=[32.6,69.6; 38,
                  69.6; 38,60; 48.8,60], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation2.y1, add3_2.u3) annotation (points=[32.6,53.6; 38,
                  53.6; 38,55.2; 48.8,55.2], style(color=74, rgbcolor={0,0,127}));
            connect(limitedAlpha.y, waveP3.u) annotation (points=[-9.5,5; -6,5;
                  -6,8; -1,8], style(color=74, rgbcolor={0,0,127}));
            connect(limitedAlpha1.y, waveP4.u) annotation (points=[-9.5,-11; -6,
                  -11; -6,-8; -1,-8], style(color=74, rgbcolor={0,0,127}));
            connect(limitedAlpha2.y, waveP5.u) annotation (points=[-9.5,-27; -6,
                  -27; -6,-24; -1,-24], style(color=74, rgbcolor={0,0,127}));
            connect(waveP3.y, fluxLimited.u1) annotation (points=[10.5,5; 14,5;
                  14,2; 19,2], style(color=74, rgbcolor={0,0,127}));
            connect(waveP4.y, fluxLimited1.u1) annotation (points=[10.5,-11; 14,
                  -11; 14,-14; 19,-14], style(color=74, rgbcolor={0,0,127}));
            connect(waveP5.y, fluxLimited2.u1) annotation (points=[10.5,-27; 14,
                  -27; 14,-30; 19,-30], style(color=74, rgbcolor={0,0,127}));
            connect(fluxLimited.y, add3_3.u1) annotation (points=[30.5,5; 34,5;
                  34,-5.4; 40.6,-5.4], style(color=74, rgbcolor={0,0,127}));
            connect(fluxLimited1.y, add3_3.u2) annotation (points=[30.5,-11; 40.6,
                  -11], style(color=74, rgbcolor={0,0,127}));
            connect(fluxLimited2.y, add3_3.u3) annotation (points=[30.5,-27; 34,
                  -27; 34,-16.6; 40.6,-16.6], style(color=74, rgbcolor={0,0,127}));
            connect(R.y, riemann.u1) annotation (points=[-79.2,57; -76,57; -76,
                  68.4; -71.2,68.4], style(color=74, rgbcolor={0,0,127}));
            connect(R.y, waveP2.u1) annotation (points=[-79.2,57; -20,57; -20,
                  52.4; -17.2,52.4], style(color=74, rgbcolor={0,0,127}));
            connect(R.y, waveP1.u1) annotation (points=[-79.2,57; -20,57; -20,
                  68.4; -17.2,68.4], style(color=74, rgbcolor={0,0,127}));
            connect(R.y, waveP.u1) annotation (points=[-79.2,57; -20,57; -20,84.4;
                  -17.2,84.4], style(color=74, rgbcolor={0,0,127}));
            connect(R.y, waveP3.u1) annotation (points=[-79.2,57; -32,57; -32,20;
                  -4,20; -4,2; -1,2], style(color=74, rgbcolor={0,0,127}));
            connect(R.y, waveP4.u1) annotation (points=[-79.2,57; -32,57; -32,20;
                  -4,20; -4,-14; -1,-14], style(color=74, rgbcolor={0,0,127}));
            connect(R.y, waveP5.u1) annotation (points=[-79.2,57; -32,57; -32,20;
                  -4,20; -4,-30; -1,-30], style(color=74, rgbcolor={0,0,127}));
            connect(lambda1.y, fluctuation.u) annotation (points=[-87.6,43; 10,43;
                  10,91.6; 18.8,91.6], style(color=74, rgbcolor={0,0,127}));
            connect(lambda2.y, fluctuation1.u) annotation (points=[-87.6,31; 10,
                  31; 10,75.6; 18.8,75.6], style(color=74, rgbcolor={0,0,127}));
            connect(lambda3.y, fluctuation2.u) annotation (points=[-87.6,19; -80,
                  19; -80,26; 12,26; 12,59.6; 18.8,59.6], style(color=74,
                  rgbcolor={0,0,127}));
            connect(lambda1.y, teta.u1) annotation (points=[-87.6,43; -74,43; -74,
                  2; -61,2],    style(color=74, rgbcolor={0,0,127}));
            connect(lambda2.y, teta1.u1) annotation (points=[-87.6,31; -76,31;
                  -76,-14; -61,-14],    style(color=74, rgbcolor={0,0,127}));
            connect(lambda3.y, teta2.u1) annotation (points=[-87.6,19; -80,19;
                  -80,-30; -61,-30],    style(color=74, rgbcolor={0,0,127}));
            connect(lambda1.y, fluxLimited.u) annotation (points=[-87.6,43; 16,43;
                  16,8; 19,8], style(color=74, rgbcolor={0,0,127}));
            connect(lambda2.y, fluxLimited1.u) annotation (points=[-87.6,31; 14,
                  31; 14,-8; 19,-8], style(color=74, rgbcolor={0,0,127}));
            connect(lambda3.y, fluxLimited2.u) annotation (points=[-87.6,19; 14,
                  19; 14,-24; 19,-24], style(color=74, rgbcolor={0,0,127}));
            connect(riemann.y, alpha.u) annotation (points=[-57.4,72; -41.2,72],
                style(color=74, rgbcolor={0,0,127}));
            connect(riemann.y, teta.u) annotation (points=[-57.4,72; -56,72; -56,
                  14; -66,14; -66,8; -61,8],    style(color=74, rgbcolor={0,0,127}));
            connect(riemann.y, teta1.u) annotation (points=[-57.4,72; -56,72; -56,
                  14; -66,14; -66,-8; -61,-8],    style(color=74, rgbcolor={0,0,
                    127}));
            connect(riemann.y, teta2.u) annotation (points=[-57.4,72; -56,72; -56,
                  14; -66,14; -66,-24; -61,-24],    style(color=74, rgbcolor={0,0,
                    127}));
            connect(alpha.y, limitedAlpha.u1) annotation (points=[-27.4,72; -26,
                  72; -26,2; -21,2], style(color=74, rgbcolor={0,0,127}));
            connect(alpha.y, limitedAlpha1.u1) annotation (points=[-27.4,72; -26,
                  72; -26,-14; -21,-14], style(color=74, rgbcolor={0,0,127}));
            connect(alpha.y, limitedAlpha2.u1) annotation (points=[-27.4,72; -26,
                  72; -26,-30; -21,-30], style(color=74, rgbcolor={0,0,127}));
            connect(deltaQ.y, riemann.u) annotation (points=[-79.6,74; -75.4,74;
                  -75.4,75.6; -71.2,75.6], style(color=74, rgbcolor={0,0,127}));
            connect(Acoustics.y, deltaQ.u)    annotation (points=[88.6,-56; 94,
                  -56; 94,98; -96,98; -96,74; -88.8,74],
                                                 style(color=74, rgbcolor={0,0,
                    127}));
            connect(add3_1.y, Acoustics.u)    annotation (points=[62.6,80; 88,
                  80; 88,-38; 28,-38; 28,-52; 31.92,-52],
                                                      style(color=74, rgbcolor={0,
                    0,127}));
            connect(add3_2.y, Acoustics.u1)    annotation (points=[62.6,60; 84,
                  60; 84,-36; 26,-36; 26,-56.8; 31.92,-56.8],
                                                      style(color=74, rgbcolor={0,
                    0,127}));
            connect(add3_3.y, Acoustics.u2)    annotation (points=[56.7,-11; 80,
                  -11; 80,-34; 22,-34; 22,-62.8; 31.92,-62.8],
                                                           style(color=74,
                  rgbcolor={0,0,127}));
            connect(teta.y, fluxLimiterSolver.u) annotation (points=[-49.5,5; -41,
                  5],                         style(color=74, rgbcolor={0,0,127}));
            connect(fluxLimiterSolver.y, limitedAlpha.u) annotation (points=[
                  -29.5,5; -24,5; -24,8; -21,8], style(color=74, rgbcolor={0,0,
                    127}));
            connect(teta1.y, fluxLimiterSolver1.u) annotation (points=[-49.5,-11;
                  -41,-11],                         style(color=74, rgbcolor={0,0,
                    127}));
            connect(fluxLimiterSolver1.y, limitedAlpha1.u) annotation (points=[
                  -29.5,-11; -24,-11; -24,-8; -21,-8], style(color=74, rgbcolor={
                    0,0,127}));
            connect(teta2.y, fluxLimiterSolver2.u) annotation (points=[-49.5,-27;
                  -41,-27],                         style(color=74, rgbcolor={0,0,
                    127}));
            connect(fluxLimiterSolver2.y, limitedAlpha2.u) annotation (points=[
                  -29.5,-27; -24,-27; -24,-24; -21,-24], style(color=74, rgbcolor=
                   {0,0,127}));
            connect(IC.y, Acoustics.u4) annotation (points=[20.4,-72; 31.92,-72],
                style(color=74, rgbcolor={0,0,127}));
            connect(BCL.y, Acoustics.u5) annotation (points=[20.4,-84; 24,-84;
                  24,-80.4; 31.92,-80.4], style(color=74, rgbcolor={0,0,127}));
            connect(BCR.y, Acoustics.u6) annotation (points=[20.4,-96; 26,-96;
                  26,-84.4; 31.92,-84.4], style(color=74, rgbcolor={0,0,127}));
          end Acoustics;

          block ICacoustics
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer m = worldModel1.m;

          equation
          for j in 1:n loop
            y[1, j] = 1.0;
            y[2, j] = 0.0;
            y[3, j] = 2.5;
          end for;

          public
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n]
              annotation (extent=[100,-10; 120,10]);
            annotation (Icon(Text(
                  extent=[-52,38; 46,-34],
                  style(color=3, rgbcolor={0,0,255}),
                  string="IC")), Documentation(info="<html>
<p>
Implements the initial condition for the acoustics example.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
          end ICacoustics;

          block BCLacoustics
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer gcl = worldModel1.gcl;

          equation
          for j in 1:gcl loop
            y[1, j] = 1.0;
            y[2, j] = 0.0;
            y[3, j] = 2.5;
          end for;

          public
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.gcl]
              annotation (extent=[100,-10; 120,10]);
            annotation (Documentation(info="<html>
<p>
Implements the left boundary condition for the acoustics example.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"), Icon(Text(
                  extent=[-48,34; 54,-38],
                  style(color=3, rgbcolor={0,0,255}),
                  string="BCL")));
          end BCLacoustics;

          block BCRacoustics
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer gcr = worldModel1.gcr;

          equation
          for j in 1:gcr loop
            y[1, j] = 0.125;
            y[2, j] = 0.0;
            y[3, j] = 0.25;
          end for;

          public
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.gcr]
              annotation (extent=[100,-10; 120,10]);
            annotation (Documentation(info="<html>
<p>
Implements the right boundary condition for the acoustics example.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"), Icon(Text(
                  extent=[-52,34; 58,-32],
                  style(color=3, rgbcolor={0,0,255}),
                  string="BCR")));
          end BCRacoustics;

          annotation (Documentation(info="<html>
<p>
This package contains acoustics example solved with the Finite Volume Methods by using flux limiters.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end Acoustic;

        package EulerSystem
          model Euler
            annotation (Diagram, Documentation(info="<html>
<h3><font color=\"#008000\" size=5>Euler equations solved with flux limiters</font></h3>
<p>
Implements the Euler system of equations
</p>
<p>
<img align=middle src=\"..\\Images\\euler.png\">
</p>

<p>
with Roe큦 flux. The initial conditions are
</p>

<img align=middle src=\"..\\Images\\sw4.png\">

<p>
The boundary conditions are
</p>

<img align=middle src=\"..\\Images\\sw6.png\">

<p>
This system has eigenvalues
</p>

<img align=middle src=\"..\\Images\\eulerroe1.png\">

<p>
and eigenvector matrix
</p>

<img align=middle src=\"..\\Images\\eulerroe2.png\">

<p>
that are computed at each time step at each interface by the corresponding blocks.
</p>


</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
            Alpha.Alpha alpha annotation (extent=[-40,66; -28,78]);
            Modelica.Blocks.Math.Add3 add3_1[worldModel1.m,worldModel1.n + 1]
              annotation (extent=[50,74; 62,86]);
            Modelica.Blocks.Math.Add3 add3_2[worldModel1.m,worldModel1.n + 1]
              annotation (extent=[50,54; 62,66]);
            Modelica.Blocks.Math.Add3 add3_3[worldModel1.m,worldModel1.n + 1]
              annotation (extent=[42,-18; 56,-4]);
            DeltaQ.deltaQ deltaQ annotation (extent=[-88,72; -80,80]);
            FLIntegrator.FLIntegrator EulerSystem(bcl=0, bcr=0)
              annotation (extent=[32,-88; 84,-48]);
            FluxLimiterSolver.FluxLimiterSolver fluxLimiterSolver
              annotation (extent=[-18,0; -8,10]);
            FluxLimiterSolver.FluxLimiterSolver fluxLimiterSolver1(p=2)
              annotation (extent=[-18,-16; -8,-6]);
            FluxLimiterSolver.FluxLimiterSolver fluxLimiterSolver2(p=3)
              annotation (extent=[-18,-32; -8,-22]);
            PDE.FiniteVolume.FluxLimiter.Examples.Acoustic.ICacoustics
              iCacoustics           annotation (extent=[12,-76; 22,-66]);
            PDE.FiniteVolume.FluxLimiter.Examples.Acoustic.BCLacoustics
              bCLacoustics            annotation (extent=[14,-84; 20,-78]);
            PDE.FiniteVolume.FluxLimiter.Examples.Acoustic.BCRacoustics
              bCRacoustics            annotation (extent=[14,-92; 20,-86]);
            Modelica.Blocks.Math.Product product[worldModel1.n + worldModel1.gcl
               + worldModel1.gcr] annotation (extent=[-86,-48; -80,-42]);
            Modelica.Blocks.Math.Division division[worldModel1.n + worldModel1.
              gcl + worldModel1.gcr] annotation (extent=[-72,-52; -66,-46]);
            Modelica.Blocks.Sources.RealExpression gamma[worldModel1.n +
              worldModel1.gcl + worldModel1.gcr](y=1.4)
              annotation (extent=[-98,-48; -92,-36]);
            Modelica.Blocks.Math.Sqrt c[worldModel1.n + worldModel1.gcl +
              worldModel1.gcr] annotation (extent=[-60,-52; -54,-46]);
            Modelica.Blocks.Math.Product product1[worldModel1.n + worldModel1.gcl
               + worldModel1.gcr] annotation (extent=[-84,-66; -78,-60]);
            Modelica.Blocks.Sources.RealExpression konst[worldModel1.n +
              worldModel1.gcl + worldModel1.gcr](y=3.5)
              annotation (extent=[-94,-64; -88,-54]);
            Modelica.Blocks.Math.Division division1[worldModel1.n + worldModel1.
              gcl + worldModel1.gcr] annotation (extent=[-72,-68; -66,-62]);
            Modelica.Blocks.Math.Product product2[worldModel1.n + worldModel1.gcl
               + worldModel1.gcr] annotation (extent=[-50,-76; -44,-70]);
            Modelica.Blocks.Math.Product product3[worldModel1.n + worldModel1.gcl
               + worldModel1.gcr] annotation (extent=[-38,-82; -32,-76]);
            Modelica.Blocks.Sources.RealExpression onehalf[worldModel1.n +
              worldModel1.gcl + worldModel1.gcr](y=0.5)
              annotation (extent=[-50,-90; -44,-78]);
            Modelica.Blocks.Math.Add h[worldModel1.n + worldModel1.gcl +
              worldModel1.gcr] annotation (extent=[-32,-68; -26,-62]);
            R r annotation (extent=[-86,28; -72,40]);
            L l(n=worldModel1.n, m=3) annotation (extent=[-88,4; -76,14]);
            WaveP1 waveP1_1 annotation (extent=[-14,82; -4,92]);
            WaveP1 waveP1_2(p=2) annotation (extent=[-14,64; -4,74]);
            WaveP1 waveP1_3(p=3) annotation (extent=[-14,46; -4,56]);
            Fluctuation1 fluctuation1_1 annotation (extent=[20,82; 30,92]);
            Fluctuation1 fluctuation1_2 annotation (extent=[20,64; 30,74]);
            Fluctuation1 fluctuation1_3 annotation (extent=[20,46; 30,56]);
            teta1 teta1_1(p=2) annotation (extent=[-36,-16; -26,-6]);
            teta1 teta1_2 annotation (extent=[-36,0; -26,10]);
            teta1 teta1_3(p=3) annotation (extent=[-36,-32; -26,-22]);
            LimitedWave1 limitedWave1_1 annotation (extent=[0,0; 10,10]);
            LimitedWave1 limitedWave1_2 annotation (extent=[0,-16; 10,-6]);
            LimitedWave1 limitedWave1_3 annotation (extent=[0,-32; 10,-22]);
            FluxLimited1 fluxLimited1_1 annotation (extent=[20,0; 30,10]);
            FluxLimited1 fluxLimited1_2(p=2) annotation (extent=[20,-16; 30,-6]);
            FluxLimited1 fluxLimited1_3(p=3) annotation (extent=[20,-32; 30,-22]);
            inner World.worldModel worldModel1(m=3,
              deltat=0.0004,
              fls=5,
              n=10)
              annotation (extent=[-100,-100; -80,-80]);
            Filter1 filter1_1 annotation (extent=[-40,38; -28,50]);
            Riemann1 riemann1_1 annotation (extent=[-68,72; -60,80]);
          equation
            connect(EulerSystem.y, deltaQ.u)  annotation (points=[86.6,-56; 94,
                  -56; 94,98; -96,98; -96,76; -88.8,76],
                                                 style(color=74, rgbcolor={0,0,
                    127}));
            connect(add3_1.y, EulerSystem.u)  annotation (points=[62.6,80; 88,
                  80; 88,-46; 28,-46; 28,-52; 29.92,-52],
                                                      style(color=74, rgbcolor={0,
                    0,127}));
            connect(add3_2.y, EulerSystem.u1)  annotation (points=[62.6,60; 84,
                  60; 84,-44; 26,-44; 26,-56.8; 29.92,-56.8],
                                                      style(color=74, rgbcolor={0,
                    0,127}));
            connect(add3_3.y, EulerSystem.u2)  annotation (points=[56.7,-11; 80,
                  -11; 80,-40; 22,-40; 22,-62.8; 29.92,-62.8],
                                                           style(color=74,
                  rgbcolor={0,0,127}));
            connect(bCRacoustics.y, EulerSystem.u6)
                                                  annotation (points=[20.3,-89;
                  25.15,-89; 25.15,-84.4; 29.92,-84.4],
                                                    style(color=74, rgbcolor={0,0,
                    127}));
            connect(bCLacoustics.y, EulerSystem.u5)
                                                  annotation (points=[20.3,-81;
                  26.15,-81; 26.15,-80.4; 29.92,-80.4],
                                                    style(color=74, rgbcolor={0,0,
                    127}));
            connect(iCacoustics.y, EulerSystem.u4)
                                                 annotation (points=[22.5,-71;
                  27.25,-71; 27.25,-72; 29.92,-72], style(color=74, rgbcolor={0,0,
                    127}));
            connect(EulerSystem.y[3, :], product.u2)
                                                   annotation (points=[86.6,-56;
                  94,-56; 94,-98; -58,-98; -58,-76; -98,-76; -98,-46.8; -86.6,
                  -46.8], style(color=74, rgbcolor={0,0,127}));
            connect(gamma.y, product.u1) annotation (points=[-91.7,-42; -88,-42;
                  -88,-43.2; -86.6,-43.2], style(color=74, rgbcolor={0,0,127}));
            connect(product.y, division.u1) annotation (points=[-79.7,-45; -75.85,
                  -45; -75.85,-47.2; -72.6,-47.2], style(color=74, rgbcolor={0,0,
                    127}));
            connect(EulerSystem.y[1, :], division.u2)
                                                    annotation (points=[86.6,-56;
                  94,-56; 94,-96; -56,-96; -56,-74; -96,-74; -96,-50.8; -72.6,
                  -50.8], style(color=74, rgbcolor={0,0,127}));
            connect(division.y, c.u) annotation (points=[-65.7,-49; -62.85,-49;
                  -62.85,-49; -60.6,-49], style(color=74, rgbcolor={0,0,127}));
            connect(konst.y, product1.u1) annotation (points=[-87.7,-59; -85.85,
                  -59; -85.85,-61.2; -84.6,-61.2], style(color=74, rgbcolor={0,0,
                    127}));
            connect(EulerSystem.y[3, :], product1.u2)
                                                    annotation (points=[86.6,-56;
                  94,-56; 94,-98; -58,-98; -58,-76; -88,-76; -88,-64.8; -84.6,
                  -64.8], style(color=74, rgbcolor={0,0,127}));
            connect(product1.y, division1.u1) annotation (points=[-77.7,-63;
                  -74.85,-63; -74.85,-63.2; -72.6,-63.2], style(color=74,
                  rgbcolor={0,0,127}));
            connect(EulerSystem.y[1, :], division1.u2)
                                                     annotation (points=[86.6,-56;
                  94,-56; 94,-96; -56,-96; -56,-74; -76,-74; -76,-66.8; -72.6,
                  -66.8], style(color=74, rgbcolor={0,0,127}));
            connect(EulerSystem.y[2, :], product2.u1)
                                                    annotation (points=[86.6,-56;
                  94,-56; 94,-94; -54,-94; -54,-71.2; -50.6,-71.2], style(color=
                    74, rgbcolor={0,0,127}));
            connect(EulerSystem.y[2, :], product2.u2)
                                                    annotation (points=[86.6,-56;
                  94,-56; 94,-94; -54,-94; -54,-74.8; -50.6,-74.8], style(color=
                    74, rgbcolor={0,0,127}));
            connect(product2.y, product3.u1) annotation (points=[-43.7,-73;
                  -40.85,-73; -40.85,-77.2; -38.6,-77.2], style(color=74,
                  rgbcolor={0,0,127}));
            connect(onehalf.y, product3.u2) annotation (points=[-43.7,-84; -42,
                  -84; -42,-80.8; -38.6,-80.8], style(color=74, rgbcolor={0,0,127}));
            connect(product3.y, h.u2) annotation (points=[-31.7,-79; -30,-79; -30,
                  -70; -36,-70; -36,-66.8; -32.6,-66.8], style(color=74, rgbcolor=
                   {0,0,127}));
            connect(division1.y, h.u1) annotation (points=[-65.7,-65; -48.85,-65;
                  -48.85,-63.2; -32.6,-63.2], style(color=74, rgbcolor={0,0,127}));
            connect(c.y, l.u1) annotation (points=[-53.7,-49; -50,-49; -50,-32; -98,-32;
                  -98,6; -89.2,6], style(color=74, rgbcolor={0,0,127}));
            connect(c.y, r.u1) annotation (points=[-53.7,-49; -50,-49; -50,-32;
                  -98,-32; -98,34; -87.4,34], style(color=74, rgbcolor={0,0,127}));
            connect(h.y, r.u2) annotation (points=[-25.7,-65; -24,-65; -24,-34;
                  -96,-34; -96,30.4; -87.4,30.4], style(color=74, rgbcolor={0,0,
                    127}));
            connect(EulerSystem.y[2, :], r.u)
                                            annotation (points=[86.6,-56; 92,-56;
                  92,-36; -94,-36; -94,37.6; -87.4,37.6], style(color=74,
                  rgbcolor={0,0,127}));
            connect(EulerSystem.y[2, :], l.u)
                                            annotation (points=[86.6,-56; 90,-56;
                  90,-36; -94,-36; -94,12; -89.2,12], style(color=74, rgbcolor={0,
                    0,127}));
            connect(alpha.y, waveP1_1.u) annotation (points=[-27.4,72; -24,72;
                  -24,90; -15,90], style(color=74, rgbcolor={0,0,127}));
            connect(alpha.y, waveP1_2.u) annotation (points=[-27.4,72; -15,72],
                style(color=74, rgbcolor={0,0,127}));
            connect(alpha.y, waveP1_3.u) annotation (points=[-27.4,72; -24,72;
                  -24,54; -15,54], style(color=74, rgbcolor={0,0,127}));
            connect(waveP1_1.y, fluctuation1_1.u1) annotation (points=[-3.5,87;
                  8.25,87; 8.25,84; 19,84], style(color=74, rgbcolor={0,0,127}));
            connect(waveP1_2.y, fluctuation1_2.u1) annotation (points=[-3.5,69;
                  7.25,69; 7.25,66; 19,66], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1_1.y, add3_1.u1) annotation (points=[30.5,89; 40,
                  89; 40,84.8; 48.8,84.8], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1_2.y, add3_1.u2) annotation (points=[30.5,71; 38,
                  71; 38,80; 48.8,80], style(color=74, rgbcolor={0,0,127}));
            connect(waveP1_3.y, fluctuation1_3.u1) annotation (points=[-3.5,51; 6,
                  51; 6,48; 19,48], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1_3.y, add3_1.u3) annotation (points=[30.5,53; 40,
                  53; 40,75.2; 48.8,75.2], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1_1.y1, add3_2.u1) annotation (points=[30.5,85; 34,
                  85; 34,64.8; 48.8,64.8], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1_2.y1, add3_2.u2) annotation (points=[30.5,67; 32,
                  67; 32,60; 48.8,60], style(color=74, rgbcolor={0,0,127}));
            connect(fluctuation1_3.y1, add3_2.u3) annotation (points=[30.5,49; 42,
                  49; 42,55.2; 48.8,55.2], style(color=74, rgbcolor={0,0,127}));
            connect(l.y, fluctuation1_1.u) annotation (points=[-75.4,12; -40,12;
                  -40,30; 10,30; 10,90; 19,90], style(color=74, rgbcolor={0,0,127}));
            connect(l.y1, fluctuation1_2.u) annotation (points=[-75.4,9; -58,9;
                  -58,16; 12,16; 12,72; 19,72], style(color=74, rgbcolor={0,0,127}));
            connect(l.y2, fluctuation1_3.u) annotation (points=[-75.4,6; -56,6;
                  -56,14; 14,14; 14,54; 19,54], style(color=74, rgbcolor={0,0,127}));
            connect(teta1_1.y, fluxLimiterSolver1.u) annotation (points=[-25.5,
                  -11; -19,-11], style(color=74, rgbcolor={0,0,127}));
            connect(l.y1, teta1_1.u1) annotation (points=[-75.4,9; -58,9; -58,-14;
                  -37,-14], style(color=74, rgbcolor={0,0,127}));
            connect(l.y, teta1_2.u1) annotation (points=[-75.4,12; -56,12; -56,2;
                  -37,2], style(color=74, rgbcolor={0,0,127}));
            connect(l.y2, teta1_3.u1) annotation (points=[-75.4,6; -56,6; -56,-30;
                  -37,-30], style(color=74, rgbcolor={0,0,127}));
            connect(teta1_3.y, fluxLimiterSolver2.u) annotation (points=[-25.5,
                  -27; -19,-27], style(color=74, rgbcolor={0,0,127}));
            connect(teta1_2.y, fluxLimiterSolver.u) annotation (points=[-25.5,5;
                  -19,5], style(color=74, rgbcolor={0,0,127}));
            connect(fluxLimiterSolver.y, limitedWave1_1.u) annotation (points=[
                  -7.5,5; -3.75,5; -3.75,8; -1,8], style(color=74, rgbcolor={0,0,
                    127}));
            connect(waveP1_1.y, limitedWave1_1.u1) annotation (points=[-3.5,87;
                  -2,87; -2,2; -1,2], style(color=74, rgbcolor={0,0,127}));
            connect(fluxLimiterSolver1.y, limitedWave1_2.u) annotation (points=[
                  -7.5,-11; -4.75,-11; -4.75,-8; -1,-8], style(color=74, rgbcolor=
                   {0,0,127}));
            connect(fluxLimiterSolver2.y, limitedWave1_3.u) annotation (points=[
                  -7.5,-27; -3.75,-27; -3.75,-24; -1,-24], style(color=74,
                  rgbcolor={0,0,127}));
            connect(waveP1_2.y, limitedWave1_2.u1) annotation (points=[-3.5,69;
                  -1,-14], style(color=74, rgbcolor={0,0,127}));
            connect(waveP1_3.y, limitedWave1_3.u1) annotation (points=[-3.5,51;
                  -1,-30], style(color=74, rgbcolor={0,0,127}));
            connect(limitedWave1_1.y, fluxLimited1_1.u1) annotation (points=[10.5,
                  5; 14.25,5; 14.25,2; 19,2], style(color=74, rgbcolor={0,0,127}));
            connect(limitedWave1_2.y, fluxLimited1_2.u1) annotation (points=[10.5,
                  -11; 14.25,-11; 14.25,-14; 19,-14], style(color=74, rgbcolor={0,
                    0,127}));
            connect(limitedWave1_3.y, fluxLimited1_3.u1) annotation (points=[10.5,
                  -27; 14.25,-27; 14.25,-30; 19,-30], style(color=74, rgbcolor={0,
                    0,127}));
            connect(l.y, fluxLimited1_1.u) annotation (points=[-75.4,12; 14,12;
                  14,8; 19,8], style(color=74, rgbcolor={0,0,127}));
            connect(l.y1, fluxLimited1_2.u) annotation (points=[-75.4,9; -68,9;
                  -68,12; 16,12; 16,-8; 19,-8], style(color=74, rgbcolor={0,0,127}));
            connect(l.y2, fluxLimited1_3.u) annotation (points=[-75.4,6; -64,6;
                  -64,12; 16,12; 16,-24; 19,-24], style(color=74, rgbcolor={0,0,
                    127}));
            connect(fluxLimited1_1.y, add3_3.u1) annotation (points=[30.5,4.9;
                  30.5,-0.55; 40.6,-0.55; 40.6,-5.4], style(color=74, rgbcolor={0,
                    0,127}));
            connect(fluxLimited1_2.y, add3_3.u2) annotation (points=[30.5,-11.1;
                  35.25,-11.1; 35.25,-11; 40.6,-11], style(color=74, rgbcolor={0,
                    0,127}));
            connect(fluxLimited1_3.y, add3_3.u3) annotation (points=[30.5,-27.1;
                  30.5,-21.55; 40.6,-21.55; 40.6,-16.6], style(color=74, rgbcolor=
                   {0,0,127}));
            connect(r.y, filter1_1.u) annotation (points=[-71.3,34; -60,34; -60,
                  44; -41.2,44], style(color=74, rgbcolor={0,0,127}));
            connect(filter1_1.y, waveP1_1.u1) annotation (points=[-27.4,44; -20,
                  44; -20,84; -15,84], style(color=74, rgbcolor={0,0,127}));
            connect(filter1_1.y, waveP1_2.u1) annotation (points=[-27.4,44; -20,
                  44; -20,66; -15,66], style(color=74, rgbcolor={0,0,127}));
            connect(filter1_1.y, waveP1_3.u1) annotation (points=[-27.4,44; -20,
                  44; -20,48; -15,48], style(color=74, rgbcolor={0,0,127}));
            connect(deltaQ.y, riemann1_1.u) annotation (points=[-79.6,76; -74,76;
                  -74,78.4; -68.8,78.4], style(color=74, rgbcolor={0,0,127}));
            connect(r.y, riemann1_1.u1) annotation (points=[-71.3,34; -68,34; -68,
                  68; -72,68; -72,73.6; -68.8,73.6], style(color=74, rgbcolor={0,
                    0,127}));
            connect(riemann1_1.y, alpha.u) annotation (points=[-59.6,76; -52,76;
                  -52,72; -41.2,72], style(color=74, rgbcolor={0,0,127}));
            connect(riemann1_1.y, teta1_2.u) annotation (points=[-59.6,76; -48,76;
                  -48,8; -37,8], style(color=74, rgbcolor={0,0,127}));
            connect(riemann1_1.y, teta1_1.u) annotation (points=[-59.6,76; -48,76;
                  -48,-8; -37,-8], style(color=74, rgbcolor={0,0,127}));
            connect(riemann1_1.y, teta1_3.u) annotation (points=[-59.6,76; -48,76;
                  -48,-24; -37,-24], style(color=74, rgbcolor={0,0,127}));
          end Euler;

          block R
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer m = worldModel1.m;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;

            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.
              gcl + worldModel1.gcr] annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + worldModel1.
              gcl + worldModel1.gcr] annotation (extent=[-140,-20; -100,20]);
            Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n + worldModel1.
              gcl + worldModel1.gcr] annotation (extent=[-140,-80; -100,-40]);
          equation
          for j in 1:n+gcl+gcr-1 loop
            y[1, 1, j] = 1.0;
            y[2, 1, j] = u[j] - u1[j];
            y[3, 1, j] = u2[j] - u[j]*u1[j];

            y[1, 2, j] = 1.0;
            y[2, 2, j] = u[j];
            y[3, 2, j] = 0.5*((u[j])^2);

            y[1, 3, j] = 1.0;
            y[2, 3, j] = u[j] + u1[j];
            y[3, 3, j] = u2[j] + u[j]*u1[j];
          end for;

          public
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.m,
              worldModel1.n + worldModel1.gcl + worldModel1.gcr - 1]
              annotation (extent=[100,-10; 120,10]);
            annotation (Diagram, Icon(
                Text(
                  extent=[-96,72; -68,52],
                  style(color=3, rgbcolor={0,0,255}),
                  string="v"),
                Text(
                  extent=[-96,12; -64,-8],
                  style(color=3, rgbcolor={0,0,255}),
                  string="c"),
                Text(
                  extent=[-94,-50; -64,-70],
                  style(color=3, rgbcolor={0,0,255}),
                  string="h")));
          end R;

          block L
            extends World.worldModel;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer m = worldModel1.m;

          equation
          for j in 1:n+1 loop
            y[j] = u[j] - u1[j];
            y1[j] = u[j];
            y2[j] = u[j] + u1[j];
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + worldModel1.gcl +
              worldModel1.gcr] annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + worldModel1.gcl +
              worldModel1.gcr] annotation (extent=[-140,-80; -100,-40]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n + 1]
              annotation (extent=[100,50; 120,70]);
            Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.n + 1]
              annotation (extent=[100,-10; 120,10]);
            Modelica.Blocks.Interfaces.RealOutput y2[worldModel1.n + 1]
              annotation (extent=[100,-70; 120,-50]);
            annotation (
              Diagram,
              Icon(Text(
                  extent=[-90,76; -66,44],
                  style(color=3, rgbcolor={0,0,255}),
                  string="v"), Text(
                  extent=[-100,-44; -60,-74],
                  style(color=3, rgbcolor={0,0,255}),
                  string="c")),
              DymolaStoredErrors);
          end L;

          block Riemann1
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer m = worldModel1.m;
          parameter Integer n = worldModel1.n;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;

          Real alpha[worldModel1.m, worldModel1.n + worldModel1.gcl + worldModel1.gcr - 1];
          equation

          for j in 1:n+gcl+gcr-1 loop
           for k in 1:m loop
                u[k, j] = u1[1, k, j]*alpha[1, j] + u1[2, k, j]*alpha[2, j] + u1[3, k, j]*alpha[3, j];
           end for;
          end for;

          y = alpha;

          // for j in 1:n+gcl+gcr-1 loop
          //   u = u1[:, :, j]*alpha;
          //   y = alpha;
          // end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
              worldModel1.gcl + worldModel1.gcr - 1]
              annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.m,
              worldModel1.n + worldModel1.gcl + worldModel1.gcr - 1]
              annotation (extent=[-140,-80; -100,-40]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + worldModel1.gcl + worldModel1.gcr - 1]
              annotation (extent=[100,-10; 120,10]);
            annotation (Diagram);
          end Riemann1;

          block WaveP1
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer m = worldModel1.m;
          parameter Integer n = worldModel1.n;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;
          parameter Integer p = 1;

          equation
          for j in 1:n+1 loop
            for i in 1:m loop
              y[i, j] = u[p, j]*u1[i, p, j];
            end for;
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
              1] annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.m,
              worldModel1.n + 1] annotation (extent=[-140,-80; -100,-40]);
            annotation (Diagram);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-10; 120,10]);
          end WaveP1;

          block Fluctuation1
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer m = worldModel1.m;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;

            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[-140,-80; -100,-40]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,30; 120,50]);
            Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-50; 120,-30]);
          equation

            for j in 1:n+1 loop
             if u[j] > 0 then
              for i in 1:m loop
                y[i, j] = u[j]*u1[i, j];
                y1[i, j] = 0;
              end for;
             elseif u[j] < 0 then
              for i in 1:m loop
                y1[i, j] = u[j]*u1[i, j];
                y[i, j] = 0;
              end for;
            else
              for i in 1:m loop
                y[i, j] = 0;
                y1[i, j] = 0;
              end for;
            end if;
           end for;

            annotation (Diagram);
          end Fluctuation1;

          block teta1
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;
          parameter Integer p = 1;

          equation
            for j in 1:n+1 loop
             if u1[j] > 0 then
              if u[p, j+1] == 0 then
                y[p, j] = 0.0;
              else
                y[p, j] = u[p, j]/u[p, j+1];
              end if;

             elseif u1[j] < 0 then
              if u[p, j+1] == 0 then
                y[p, j] = 0.0;
              else
                y[p, j] = u[p, j+2]/u[p, j+1];
              end if;

             else
               y[p, j] = 0.0;
             end if;
            end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
              worldModel1.gcl + worldModel1.gcr - 1]
              annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
              annotation (extent=[-140,-80; -100,-40]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-10; 120,10]);
            annotation (Diagram);
          end teta1;

          block LimitedWave1
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer m = worldModel1.m;
          parameter Integer n = worldModel1.n;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;
          parameter Integer p = 1;

          equation
          for j in 1:n+1 loop
            for i in 1:m loop
              y[i, j] = u[p, j]*u1[i, p];
            end for;
          end for;

          public
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-10; 120,10]);
            annotation (Diagram);
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.n +
              1] annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[-140,-80; -100,-40]);
          end LimitedWave1;

          block FluxLimited1
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer n = worldModel1.n;
          parameter Integer m = worldModel1.m;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;
          parameter Real deltat = worldModel1.deltat;
          parameter Real deltax = 1/n;
          parameter Integer p = 1;

          equation
          for j in 1:n+1 loop
            for i in 1:m loop
              y[i, j] = 0.5*abs(u[j])*(1-(deltat/deltax)*abs(u[j]))*u1[i, j];
            end for;
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.n + 1]
              annotation (extent=[-140,40; -100,80]);
            Modelica.Blocks.Interfaces.RealInput u1[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[-140,-80; -100,-40]);
            annotation (Diagram);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.n
               + 1] annotation (extent=[100,-12; 120,8]);
          end FluxLimited1;

          block Filter1
            extends Icons.BlockIcon;

          outer PDE.World.worldModel worldModel1;
          parameter Integer m = worldModel1.m;
          parameter Integer n = worldModel1.n;
          parameter Integer gcl = worldModel1.gcl;
          parameter Integer gcr = worldModel1.gcr;
          parameter Integer p = 1;

          equation
          for j in 1:n+1 loop
            for i in 1:m loop
              for k in 1:m loop
                y[i, k, j] = u[i, k, j+1];
              end for;
            end for;
          end for;

          public
            Modelica.Blocks.Interfaces.RealInput u[worldModel1.m,worldModel1.m,
              worldModel1.n + worldModel1.gcl + worldModel1.gcr - 1]
              annotation (extent=[-140,-20; -100,20]);
            Modelica.Blocks.Interfaces.RealOutput y[worldModel1.m,worldModel1.m,
              worldModel1.n + 1] annotation (extent=[100,-10; 120,10]);
            annotation (Diagram);
          end Filter1;
          annotation (Documentation(info="<html>
<p>
This package contains Euler equations example solved with the Finite Volume Methods by using flux limiters.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        end EulerSystem;
        annotation (Documentation(info="<html>
<p>
This package contains examples of partial differential equations solved with the Finite Volume Methods by using flux limiters.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end Examples;
      annotation (Documentation(info="<html>
<p>
The Flux Limiter package offers all the necessary components to solve a linear constant-coefficient system of hyperbolic equations. <br>
In this package <b>m</b> is the number of equations in the system, <b>n</b> is the number of unknown cell averages, <b>gcl</b> and <b>gcr</b> <br>
are the ghost cells at the left and at the right part of the domain respectively. <br>
For the explanations of how to use the blocks in this package see the <b>Users Guide</b>.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end FluxLimiter;

    annotation (Documentation(info="<html>
<p>
This package contains necessary blocks for solving partial differential equations with Finite Volume Methods.
To understand the use of the blocks, many examples are implemented (PDE->FiniteVolume->Examples).
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
  end FiniteVolume;

  package General
    package Integrator
      block GeneralIntegrator
        extends Icons.BlockIcon1;

        outer PDE.World.worldModel worldModel1;
        inner parameter Integer n = worldModel1.n;
        parameter Integer integrator = 1;

        inner parameter Integer bcl = 0
          "|Boundary Conditions| Type of the boundary condition at the left (-1:symmtery; 0: none)";
        inner parameter Integer bcr = 0
          "|Boundary Conditions| Type of the boundary condition at the right (-1:symmtery; 0: none)";
        inner parameter Integer vb = 1;
        inner parameter Integer ve = n;
        inner parameter Integer icb = 1;
        inner parameter Integer ice = n;

        Modelica.Blocks.Interfaces.RealInput u[worldModel1.n]
          annotation (extent=[-118,70; -100,94]);
        Modelica.Blocks.Interfaces.RealInput u1[worldModel1.n + 1]
          annotation (extent=[-118,28; -100,52]);
        Modelica.Blocks.Interfaces.RealInput u2[worldModel1.n]
          annotation (extent=[-118,-30; -100,-6]);
        Modelica.Blocks.Interfaces.RealInput u3
          annotation (extent=[-118,-60; -100,-36]);
        Modelica.Blocks.Interfaces.RealInput u4
          annotation (extent=[-118,-90; -100,-66]);
        annotation (Diagram, Icon(
            Text(
              extent=[-92,92; -44,70],
              style(color=3, rgbcolor={0,0,255}),
              string="R(MOL)"),
            Text(
              extent=[-92,48; -46,30],
              style(color=3, rgbcolor={0,0,255}),
              string="F(FVM)"),
            Text(
              extent=[-106,-12; -58,-28],
              style(color=3, rgbcolor={0,0,255}),
              string="IC"),
            Text(
              extent=[-96,-42; -58,-58],
              style(color=3, rgbcolor={0,0,255}),
              string="BCL"),
            Text(
              extent=[-94,-72; -56,-88],
              style(color=3, rgbcolor={0,0,255}),
              string="BCR"),
            Text(
              extent=[52,70; 94,50],
              style(color=3, rgbcolor={0,0,255}),
              string="Var"),
            Text(
              extent=[62,-50; 90,-70],
              style(color=3, rgbcolor={0,0,255}),
              string="Q")),
          Documentation(info="<html>
<p>
The GeneralIntegrator block combines two integrators into one block. The scope of the block is to make more abstract the <br>
procedure of solving PDE. This provide transparence to the user, that can just take one integrator block and choose the numerical <br>
method to use. The integrator parameter allows to select the method: 1 for MOL and 2 for FVM. The treatment of the IC, BCL and BCR <br>
inputs is the same as for other integrator blocks.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        Modelica.Blocks.Interfaces.RealOutput y[worldModel1.n]
          annotation (extent=[100,50; 120,70]);
        Modelica.Blocks.Interfaces.RealOutput y1[worldModel1.n + worldModel1.gcl +
          worldModel1.gcr] annotation (extent=[100,-70; 120,-50]);
        PDE.General.Interfaces.IntegratorMOL integratorMOL if
                                       integrator == 1
          annotation (extent=[-18,28; 28,76]);
        PDE.General.Interfaces.IntegratorFVM integratorFVM(
                                    bcl=0, bcr=0) if
                                       integrator == 2
          annotation (extent=[-18,-60; 28,-18]);
      equation
        connect(u, integratorMOL.u) annotation (points=[-109,82; -70,82; -70,
              66.88; -20.07,66.88], style(color=74, rgbcolor={0,0,127}));
        connect(u1, integratorFVM.u) annotation (points=[-109,40; -70,40; -70,
              -22.2; -20.53,-22.2], style(color=74, rgbcolor={0,0,127}));
        connect(u2, integratorMOL.u1) annotation (points=[-109,-18; -80,-18; -80,
              49.6; -20.07,49.6], style(color=74, rgbcolor={0,0,127}));
        connect(u2, integratorFVM.u2) annotation (points=[-109,-18; -80,-18; -80,
              -39; -20.53,-39], style(color=74, rgbcolor={0,0,127}));
        connect(u3, integratorMOL.u2) annotation (points=[-109,-48; -50,-48; -50,
              42.4; -20.07,42.4], style(color=74, rgbcolor={0,0,127}));
        connect(u3, integratorFVM.u3[1]) annotation (points=[-109,-48; -70.265,
              -48; -70.265,-50.97; -20.53,-50.97], style(color=74, rgbcolor={0,0,
                127}));
        connect(u3, integratorFVM.u3[2]) annotation (points=[-109,-48; -70.265,
              -48; -70.265,-48.03; -20.53,-48.03], style(color=74, rgbcolor={0,0,
                127}));
        connect(u4, integratorFVM.u4[1]) annotation (points=[-109,-78; -40,-78;
              -40,-57.27; -20.53,-57.27], style(color=74, rgbcolor={0,0,127}));
        connect(u4, integratorFVM.u4[2]) annotation (points=[-109,-78; -40,-78;
              -40,-54.33; -20.53,-54.33], style(color=74, rgbcolor={0,0,127}));
        connect(u4, integratorMOL.u3) annotation (points=[-109,-78; -40,-78; -40,
              35.2; -20.07,35.2], style(color=74, rgbcolor={0,0,127}));
        connect(integratorMOL.y, y) annotation (points=[30.3,66.4; 66.15,66.4;
              66.15,60; 110,60], style(color=74, rgbcolor={0,0,127}));
        connect(integratorFVM.y, y1) annotation (points=[30.3,-26.4; 60,-26.4; 60,
              -60; 110,-60], style(color=74, rgbcolor={0,0,127}));
      end GeneralIntegrator;
      annotation (Documentation(info="<html>
<p>
This package contains GeneralIntegrator block that implements the Method of Lines and Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end Integrator;

    package Interfaces
      block IntegratorMOL
        extends PDE.Icons.BlockIcon1;

      outer Integer n;
      outer Integer bcl;
      outer Integer bcr;

      outer Integer vb;
      outer Integer ve;
      outer Integer icb;
      outer Integer ice;

      Real f[n];

      equation
        y = f;

        if bcl == 1 then
          f[1] = u2;
        end if;
        if bcr == 1 then
          f[n] = u3;
        end if;

        for i in vb:ve loop
          der(f[i]) = u[i];
        end for;

      initial equation

        for i in icb:ice loop
          f[i] = u1[i];
        end for;

        annotation (Icon(
            Text(
              extent=[-14,126; 16,80],
              style(color=3, rgbcolor={0,0,255}),
              string="%name"),
            Text(
              extent=[-98,-2; -68,-20],
              style(color=3, rgbcolor={0,0,255}),
              string="IC"),
            Text(
              extent=[-96,-32; -60,-48],
              style(color=3, rgbcolor={0,0,255}),
              string="BCL"),
            Text(
              extent=[-98,-62; -54,-78],
              style(color=3, rgbcolor={0,0,255}),
              string="BCR"),
            Text(
              extent=[-96,72; -68,48],
              style(color=3, rgbcolor={0,0,255}),
              string="R"),
            Text(
              extent=[54,74; 94,48],
              style(color=3, rgbcolor={0,0,255}),
              string="Var")), Diagram,
          Documentation(info="<html>
<p>
This is a slightly modified version of the MOL integrator block implemented in PDE->MOL->Integrator. The integrator block accepts the equations of the form
</p>
<p>
u_t = R(u, u_xx, ...)
</p>
<p>
where R is a function of u, u_xx, and so on. Once the equation is transformed in this form <br>
we can construct the right part R and pass it to the input R of the integrator block. <br>
The unknown varible u which we need for the construction of the right part of the equation is <br>
provided by the Var output of the integrator block.
The initial condition is passed to the IC input, whereas the left and right boundary conditions are <br>
passed to the BCL and BCR inputs respectively.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      public
        Modelica.Blocks.Interfaces.RealInput u[n]
          annotation (extent=[-118,50; -100,74]);
        Modelica.Blocks.Interfaces.RealOutput y[n]
          annotation (extent=[100,50; 120,70]);
        Modelica.Blocks.Interfaces.RealInput u1[n]
          annotation (extent=[-118,-22; -100,2]);
        Modelica.Blocks.Interfaces.RealInput u2
          annotation (extent=[-118,-52; -100,-28]);
        Modelica.Blocks.Interfaces.RealInput u3
          annotation (extent=[-118,-82; -100,-58]);
      equation

      end IntegratorMOL;

      block IntegratorFVM

        extends Icons.BlockIcon1;

      outer Integer n;

      parameter Integer vb = gcl+1 "|Unknowns| The left most unknown";
      parameter Integer ve = gcl + n "|Unknowns| The right most unknown";

      parameter Integer icb = gcl+1
          "|Initial Condition| Begin of the initial condition";
      parameter Integer ice = gcl + n
          "|Initial Condition| End of the initial condition";

      parameter Integer bcl = 1
          "|Boundary Conditions| Boundary condition at the left (0: no; 1: yes)";
      parameter Integer bcr = 1
          "|Boundary Conditions| Boundary condition at the right (0: no; 1: yes)";

      parameter Integer gcl = 2
          "|Boundary Conditions| Number of ghost cells at the left";
      parameter Integer gcr = 2
          "|Boundary Conditions| Number of ghost cells at the right";

      parameter Real delta_x = 1/n;
      //parameter Real delta_t = 0.1;

      Real q[n+gcl+gcr];

      equation
        y = q;

        for i in 1:gcl loop
           q[i]= u3[i];
        end for;

        if bcl == 1 then
           q[gcl+1] = q[gcl];
        end if;

        for i in 1:gcr loop
           q[gcl+n+i] = u4[i];
        end for;

        if bcr == 1 then
           q[gcl+n] = q[gcl+n+1];
        end if;

        for i in vb:ve loop
          der(q[i]) = -(1/delta_x)*(u[i-gcl+1]-u[i-gcl]);
        end for;

      initial equation

        for i in icb:ice loop
          q[i] =  u2[i-gcl];
        end for;

      // equation
      //   y = q;
      //
      //   if bcl == 1 then
      //     for i in 1:gcl loop
      //       q[i]= u3[i];
      //     end for;
      //   end if;
      //
      //   if bcr == 1 then
      //     for i in 1:gcr loop
      //       q[gcl+n+i] = u4[i];
      //     end for;
      //   end if;
      //
      //   for i in vb:ve loop
      //     der(q[i]) = -(1/delta_x)*(u[i-gcl+1]-u[i-gcl]);
      //   end for;
      //
      // initial equation
      //
      //   for i in icb:ice loop
      //     q[i] =  u2[i-gcl];
      //   end for;

      public
        Modelica.Blocks.Interfaces.RealInput u[n + 1]
          annotation (extent=[-122,66; -100,94]);
        Modelica.Blocks.Interfaces.RealOutput y[n + 4]
          annotation (extent=[100,50; 120,70]);
        Modelica.Blocks.Interfaces.RealInput u2[n]
          annotation (extent=[-122,-14; -100,14]);
        Modelica.Blocks.Interfaces.RealInput u3[2]
          annotation (extent=[-122,-64; -100,-36]);
        Modelica.Blocks.Interfaces.RealInput u4[2]
          annotation (extent=[-122,-94; -100,-66]);
        annotation (Diagram, Icon(
            Text(
              extent=[-40,162; 40,124],
              style(color=3, rgbcolor={0,0,255}),
              string="%name"),
            Text(
              extent=[-108,88; -60,70],
              style(color=3, rgbcolor={0,0,255}),
              string="F"),
            Text(
              extent=[-100,10; -64,-10],
              style(color=3, rgbcolor={0,0,255}),
              string="IC"),
            Text(
              extent=[-102,-38; -58,-58],
              style(color=3, rgbcolor={0,0,255}),
              string="gcl"),
            Text(
              extent=[-108,-68; -50,-88],
              style(color=3, rgbcolor={0,0,255}),
              string="gcr"),
            Text(
              extent=[60,74; 90,46],
              style(color=3, rgbcolor={0,0,255}),
              string="Q")),
          Documentation(info="<html>
<p>
The IntegratorFVM block is a slightly modified version of the FVMIntegrator (PDE->FiniteVolume->FVMIntegrator).
It implements the cell average update rule. In one dimension, the finite volume method subdivide the domain into <br>
cells (intervals) and approximates the integral of the unknown function q over each of these cells at each time step (see figure below). <br>
The ghost cells are the boundary cells that are introduced to avoid writing special formulas for the boundary cells.
</p>

<img align=middle src=\"..\\Images\\fvmi7.png\">

<p>
Denote the i-th cell by
</p>

<img align=middle src=\"..\\Images\\fvmi.png\">

<p>
Then the approximation to the average of q in the cell C<sub>i</sub> at time t, which we denote with Q<sub>i</sub> is
</p>

<img align=middle src=\"..\\Images\\fvmi1.png\">


<p>
The approximation to this average derives from the integral form of the conservation law
</p>

<img align=middle src=\"..\\Images\\fvmi2.png\">

<p>
which states that the average within the cell can only changes due to the fluxes at the boundaries (if we assume that <br>
no source or sink is present in the cell). <br>
If we integrate this expression in time from t to t+deltat, we obtain <br>
</p>

<img align=middle src=\"..\\Images\\fvmi3.png\">

<p>
and dividing by deltax we reach the form
</p>

<img align=middle src=\"..\\Images\\fvmi4.png\">

<p>
which give us an explicit time marching algorithm. By using the notation for averages introduced above we can write
</p>

<img align=middle src=\"..\\Images\\fvmi5.png\">

<p>
where
</p>

<img align=middle src=\"..\\Images\\fvmi6.png\">

<p>
approximates the average flux along the interface x<sub>i-1/2</sub>. <br>
The FVMIntegrator block implements this average update rule. Initial condition of the problem can be passed to the IC input of the block, <br>
whereas the boundary conditions to the corresponding ghost cells. The number of ghost cells depends on the method we use. In the present <br>
package, two ghost cells at the left and at the right are enough for every method implemented.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
      end IntegratorFVM;
      annotation (Documentation(info="<html>
<p>
This package contains two integrator blocks: IntegratorMOL and IntegratorFVM.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end Interfaces;

    package Examples
      model AdvectionMOL
        PDE.General.Integrator.GeneralIntegrator generalIntegrator(
          vb=2,
          icb=2,
          n=worldModel1.n,
          ve=worldModel1.n,
          ice=worldModel1.n,
          bcl=1,
          integrator=1) annotation (extent=[8,-16; 50,24]);
        MOL.SpaceDerivative.Derivatives.u_x u_x
          annotation (extent=[-58,20; -38,40]);
        Modelica.Blocks.Math.Product product[worldModel1.n]
          annotation (extent=[-26,14; -14,26]);
        annotation (Diagram, Documentation(info="<html>
<p>
Implements advection equation with MOL method by using GeneralIntegrator block and setting the integrator parameter to 1.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
        Modelica.Blocks.Sources.RealExpression speed[worldModel1.n](y=-0.1)
          annotation (extent=[-52,2; -44,18]);
        MOL.Examples.Diffusion.DiffusionIC diffusionIC
          annotation (extent=[-26,-6; -14,6]);
        Modelica.Blocks.Sources.RealExpression BCL(y=cos(-0.1*time))
          annotation (extent=[-30,-28; -10,-8]);
        inner World.worldModel worldModel1 annotation (extent=[-20,88; 20,100]);
      equation
        connect(u_x.y, product.u1) annotation (points=[-36.9,30; -32,30; -32,
              23.6; -27.2,23.6],
                           style(color=74, rgbcolor={0,0,127}));
        connect(speed.y, product.u2) annotation (points=[-43.6,10; -34,10; -34,
              16.4; -27.2,16.4], style(color=74, rgbcolor={0,0,127}));
        connect(product.y, generalIntegrator.u) annotation (points=[-13.4,20;
              -3.645,20; -3.645,20.4; 6.11,20.4],
                         style(color=74, rgbcolor={0,0,127}));
        connect(generalIntegrator.y, u_x.u) annotation (points=[52.1,16; 60,16;
              60,48; -66,48; -66,30; -60,30], style(color=74, rgbcolor={0,0,127}));
        connect(diffusionIC.y, generalIntegrator.u2) annotation (points=[-13.4,0;
              -3.645,0; -3.645,0.4; 6.11,0.4], style(color=74, rgbcolor={0,0,127}));
        connect(BCL.y, generalIntegrator.u3) annotation (points=[-9,-18; -4,-18;
              -4,-5.6; 6.11,-5.6],
                                 style(color=74, rgbcolor={0,0,127}));
      end AdvectionMOL;

      model AdvectionFVM
        PDE.General.Integrator.GeneralIntegrator generalIntegrator(
          vb=2,
          icb=2,
          n=worldModel1.n,
          ve=worldModel1.n,
          ice=worldModel1.n,
          bcl=0,
          integrator=2) annotation (extent=[-28,16; 14,56]);
        annotation (Diagram, Documentation(info="<html>
<p>
Implements advection equation with FVM method by using GeneralIntegrator block and setting the integrator parameter to 2.
</p>
</html>"));
        MOL.Examples.Diffusion.DiffusionIC diffusionIC
          annotation (extent=[-68,28; -56,40]);
        Modelica.Blocks.Sources.RealExpression BCL(y=cos(-0.1*time))
          annotation (extent=[-72,6; -52,26]);
        inner World.worldModel worldModel1 annotation (extent=[-28,88; 14,100]);
        PDE.FiniteVolume.LDLR.L.LDLRminus lDLRminus_plus
          annotation (extent=[32,28; 54,36]);
        PDE.FiniteVolume.LDLR.L.LDLRplus lDLRplus_plus
          annotation (extent=[32,16; 54,24]);
        PDE.FiniteVolume.LDLR.u.u_minus u_minus_plus
          annotation (extent=[66,28; 82,36]);
        PDE.FiniteVolume.LDLR.u.u_plus u_plus_plus
          annotation (extent=[66,16; 82,24]);
        Modelica.Blocks.Sources.RealExpression speed1[worldModel1.n + 1](y=0.1)
          annotation (extent=[10,-24; 24,-8]);
        Modelica.Blocks.Math.Product product1[worldModel1.n + 1]
          annotation (extent=[38,-12; 50,0]);
        Modelica.Blocks.Math.Product product2[worldModel1.n + 1]
          annotation (extent=[38,-32; 50,-20]);
        PDE.FiniteVolume.Fluxes.LaxFriedrichFlux.LF lF(
                                        alpha=0.2)
          annotation (extent=[64,-26; 84,-6]);
        Modelica.Blocks.Sources.RealExpression BCR
          annotation (extent=[-68,-10; -56,10]);
      equation
        connect(diffusionIC.y, generalIntegrator.u2) annotation (points=[-55.4,34;
              -41.75,34; -41.75,32.4; -29.89,32.4],
                                               style(color=74, rgbcolor={0,0,127}));
        connect(BCL.y, generalIntegrator.u3) annotation (points=[-51,16; -40,16;
              -40,26.4; -29.89,26.4],
                                 style(color=74, rgbcolor={0,0,127}));
        connect(generalIntegrator.y1, lDLRminus_plus.u) annotation (points=[16.1,24;
              24,24; 24,32; 29.8,32],     style(color=74, rgbcolor={0,0,127}));
        connect(generalIntegrator.y1, lDLRplus_plus.u) annotation (points=[16.1,24;
              24,24; 24,20; 29.8,20],     style(color=74, rgbcolor={0,0,127}));
        connect(lDLRminus_plus.y, u_minus_plus.u1) annotation (points=[55.1,32;
              60,32; 60,29.6; 64.4,29.6], style(color=74, rgbcolor={0,0,127}));
        connect(lDLRplus_plus.y, u_plus_plus.u1) annotation (points=[55.1,20; 60,
              20; 60,17.6; 64.4,17.6], style(color=74, rgbcolor={0,0,127}));
        connect(generalIntegrator.y1, u_minus_plus.u) annotation (points=[16.1,24;
              24,24; 24,40; 60,40; 60,34.4; 64.4,34.4], style(color=74, rgbcolor=
                {0,0,127}));
        connect(generalIntegrator.y1, u_plus_plus.u) annotation (points=[16.1,24;
              24,24; 24,26; 60,26; 60,22.4; 64.4,22.4], style(color=74, rgbcolor=
                {0,0,127}));
        connect(u_plus_plus.y, product2.u2) annotation (points=[82.8,20; 86,20;
              86,8; 4,8; 4,-29.6; 36.8,-29.6],   style(color=74, rgbcolor={0,0,
                127}));
        connect(u_minus_plus.y, product1.u1) annotation (points=[82.8,32; 88,32;
              88,6; 30,6; 30,-2.4; 36.8,-2.4], style(color=74, rgbcolor={0,0,127}));
        connect(speed1.y, product1.u2) annotation (points=[24.7,-16; 30,-16; 30,
              -9.6; 36.8,-9.6], style(color=74, rgbcolor={0,0,127}));
        connect(speed1.y, product2.u1) annotation (points=[24.7,-16; 30,-16; 30,
              -22.4; 36.8,-22.4], style(color=74, rgbcolor={0,0,127}));
        connect(product1.y, lF.u) annotation (points=[50.6,-6; 56,-6; 56,-8; 63.1,
              -8], style(color=74, rgbcolor={0,0,127}));
        connect(product2.y, lF.u1) annotation (points=[50.6,-26; 54,-26; 54,-11;
              63.1,-11], style(color=74, rgbcolor={0,0,127}));
        connect(u_plus_plus.y, lF.u2) annotation (points=[82.8,20; 86,20; 86,8;
              60,8; 60,-21; 63.1,-21], style(color=74, rgbcolor={0,0,127}));
        connect(u_minus_plus.y, lF.u3) annotation (points=[82.8,32; 88,32; 88,6;
              58,6; 58,-24; 63.1,-24], style(color=74, rgbcolor={0,0,127}));
        connect(BCR.y, generalIntegrator.u4) annotation (points=[-55.4,0; -38,0;
              -38,20.4; -29.89,20.4],
                                 style(color=74, rgbcolor={0,0,127}));
        connect(lF.y, generalIntegrator.u1) annotation (points=[85,-16; 90,-16;
              90,-40; -86,-40; -86,44; -29.89,44], style(color=74, rgbcolor={0,0,
                127}));
      end AdvectionFVM;
      annotation (Documentation(info="<html>
<p>
This package contains examples of partial differential equations solved with the general block.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
    end Examples;
    annotation (Documentation(info="<html>
<p>
This package contains necessary blocks for solving partial differential equations with two numerical methods: Method of Lines and Finite Volume Methods.
</p>

</pre>
<p><b>Release Notes: </b></p>

<ul>
</html>"));
  end General;
end PDE;
