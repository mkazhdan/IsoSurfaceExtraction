<CENTER><H1>Accurate Isosurface Interpolation with Hermite Data</A></H1></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#EXECUTABLE">executable</A>
</CENTER>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A href="http://www.cs.jhu.edu/~misha/MyPapers/3DV15.pdf">3DV 2015 Paper</A><br>
<A href="IsoSurfaceExtraction.x64.exe">Windows Executable (x64)</A><br>
<A href="IsoSurfaceExtraction.zip">Source Code</A><br>
<A href="license.txt">License</A><br>
<HR>
<A NAME="DESCRIPTION"><B>CODE DESCRIPTION</B></A><br>
<UL>
The code provides an implementation of <A HREF="http://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Lorensen87.pdf">Marching-Cubes</A> isosurface extraction, supporting the use of Hermite Data and a full case table for resolving the ambiguous cases in the original implementation of Lorensen and Cline. Given a regular grid of values and a prescribed isovalue, the corresponding isosurface is extracted and output to file.
<P>
Implementations of Hermite interpolation on adapted octrees can also be found in the code for <A HREF="http://www.cs.jhu.edu/~misha/Code/PoissonRecon">Poisson Surface Reconstruction</A>, <A HREF="http://www.gris.informatik.tu-darmstadt.de/projects/floating-scale-surface-recon/">Floating Scale Surface Reconstruction</A>, and <A HREF="http://www.cs.jhu.edu/~misha/Code/IsoOctree/">Unconstrained Isosurface Extraction on Arbitrary Octrees</A>.
</UL>
<HR>
<A NAME="EXECUTABLE"><B>EXECUTABLE ARGUMENTS</B></A><br>
<UL>
<DL>

<DT><b>--in</b> &#60;<i>input volume</i>&#62;
<DD> This string is the name of the file containing the regularly sampled voxeld values. Values are linearized, so that if the image is of resolution <I>ResX</I>*<I>ResY</I>*<I>ResY</I> then the value at index (<I>x,y,z</I>) is at position <I>x</I>+<I>y</I>*<I>ResX</I>+<I>z</I>*<I>ResX</I>*<I>ResY</I>, and are stored in binary format either as unsigned chars (1 byte per value)  or as floating point values (4 bytes per value).

<DT><b>--res</b> &#60;<i>x-resolution y-resolution z-resolution</i>&#62;
<DD> This triplet of integer values gives the resolution of the voxel grid.

<DT>[<b>--dim</b> &#60;<i>x-scale y-scale z-scale</i>&#62;]
<DD> This (white-space separated) triplet of floating point values gives the scale of a voxel. (By default, the values are set to 1.)

<DT>[<b>--out</b> &#60;<i>output mesh</i>&#62;]
<DD> This string is the name of the file to which the extracted isosurface will be written. The mesh will be written in the
<A HREF="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--iso</b> &#60;<i>isovalue</i>&#62;]
<DD> This floating point value specifies the isovalue at which the isosurface is to be extracted. (By default, the value is set to 0.)
  
<DT>[<b>--sIters</b> &#60;<i>smoothing iterations</i>&#62;]
<DD> This integer specified the number of (1-ring) smoothing iterations that are to be applied to the voxel grid before extracting the isosurface. (By default, the value is set to 0.)

<DT>[<b>--full</b>]
<DD> If specified, the Marching-Cubes algorithm is implemented using a ``full'' case table, using the average value of face corners to resolve the amiguous case when the face is zero-crossing and the values on antipodal corners are the same.

<DT>[<b>--flip</b>]
<DD> If specified, the triangle orientation in the output is flipped.

<DT>[<b>--quadratic</b>]
<DD> If specified, Hermite data is generated (by considering the difference between neighboring values) and a second-order interpolant is used to define the positions of zero-crossings.

<DT>[<b>--polygon</b>]
<DD> If specified, the raw iso-polygons are output. Otherwise, the polygons are triangulated using a minimal-area-triangulation and a triangle mesh is written as output.
    
<DT>[<b>--float</b>]
<DD> If specified, the input is read in as a list of floating point values. Otherwise, the input is assumed to represent (unsigned) char values.
  
<DT>[<b>--nonManifold</b>]
<DD> Although the polygon mesh resulting from Marching-Cubes is manifold, it can be the case that two polygon share two vertices that are not on a polygon edge. As a result, the minimal area triangulation could introduce the same edge for the triangulation of both polygons, resulting in a triangle mesh with non-manifold edges. By default, the code will introduce an additional vertex (the plane's barycenter) if a non-manifold triangulation can arise, thereby ensuring that the output triangle mesh is manifold. Enabling this flag avoids introducing the barycenter, but could result in mesh with non-manifold edges.

</UL>
<HR>
<A HREF="http://www.cs.jhu.edu/~misha">HOME</A>
