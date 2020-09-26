using System.Numerics;
using UnityEngine;

public class streamlines : MonoBehaviour
{
    public solver solver;
    public MeshFilter streamline0;
    public MeshFilter streamline1;
    public MeshFilter streamline2;
    public MeshFilter streamline3;
    public MeshFilter streamline4;
    public MeshFilter streamline5;
    public MeshFilter streamline6;
    public MeshFilter streamline7;
    public MeshFilter streamline8;
    public MeshFilter streamline9;
    float xDisplacement = 0;
    float yDisplacement = 0;
    float angleOfAttack = 0;

    UnityEngine.Vector3[] vertices, nv, newVertices, zeroVertices;

    int[] indicesForStreamlines = new int[101];
    int[] zeroIndices = new int[101];
    [HideInInspector]public enum DisplayGrid { joukowsky, nacaFourDigit, nacaFourDigitModified, nacaFiveDigit, nacaFiveDigitModified };
    int airfoilType = -1;

    Complex J(Complex x, float y, float z)
    {
        return x + (Complex.Exp(-Complex.ImaginaryOne * 2 * z) * y * y) / x;
    }

    void Joukowsky()
    {
        if (solver.airfoilType == solver.DisplayGrid.joukowsky && solver.displayJoukowsky != solver.DisplayJoukowsky.numericSolver)
        {
            if (solver.xDisplacement != xDisplacement || solver.yDisplacement != yDisplacement || solver.angleOfAttack != angleOfAttack || (int)solver.airfoilType != airfoilType)
            {
                xDisplacement = solver.xDisplacement;
                yDisplacement = solver.yDisplacement;
                angleOfAttack = solver.angleOfAttack;
                airfoilType = (int)solver.airfoilType;
                float alpha = solver.angleOfAttack * Mathf.PI / 180f * 0.035f;
                float k = 2f * solver.freestreamVelocity * Mathf.Sin(alpha);
                float circulation = k / (2f * Mathf.PI);
                Complex s = new Complex(-solver.xDisplacement, solver.yDisplacement);
                float lambda = 1 + (float)s.Magnitude;
                Complex flow, w, zeta;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == -10)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.01 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary + 7.5f);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline0.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline0.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline0.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline0.mesh.vertices = newVertices;
                streamline0.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline0.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == -8)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.02 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary + 6);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline1.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline1.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline1.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline1.mesh.vertices = newVertices;
                streamline1.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline1.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == -6)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.03 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary + 4.5f);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline2.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline2.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline2.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline2.mesh.vertices = newVertices;
                streamline2.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline2.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == -4)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.06 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary + 3);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline3.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline3.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline3.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline3.mesh.vertices = newVertices;
                streamline3.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline3.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == -2)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.2 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary + 1.5f);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline4.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline4.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline4.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline4.mesh.vertices = newVertices;
                streamline4.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline4.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == 2)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.2 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary - 1.5f);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline5.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline5.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline5.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline5.mesh.vertices = newVertices;
                streamline5.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline5.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == 4)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.06 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary - 3);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline6.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline6.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline6.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline6.mesh.vertices = newVertices;
                streamline6.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline6.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == 6)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.03 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary - 4.5f);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline7.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline7.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline7.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline7.mesh.vertices = newVertices;
                streamline7.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline7.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == 8)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.02 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary - 6);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline8.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline8.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline8.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline8.mesh.vertices = newVertices;
                streamline8.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline8.GetComponent<Renderer>().enabled = true;

                for (int i = -10; i <= 10; ++i)
                {
                    for (int j = -50; j <= 50; ++j)
                    {
                        if (i == 10)
                        {
                            nv = new UnityEngine.Vector3[101];
                            Complex z = new Complex(j, i);
                            w = J(z, lambda, alpha);
                            if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                            else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                            else if (j == 0) zeta = new Complex(j, 1.01 * i);
                            flow = -zeta * Complex.Exp(Complex.ImaginaryOne * alpha) + (Complex.Exp(-Complex.ImaginaryOne * alpha) * Complex.Abs(s - 1) * Complex.Abs(s - 1) / (zeta - s)) + (Complex.ImaginaryOne * k * Complex.Log(zeta - s));
                            newVertices[j + 50] = new UnityEngine.Vector3(2.5f * (float)w.Real, -2.1f * (float)flow.Imaginary - 7.5f);
                            indicesForStreamlines[j + 50] = j + 50;
                        }
                    }
                }

                streamline9.transform.localScale = new UnityEngine.Vector3(solver.chord * 0.1f, solver.chord * 0.1f, solver.chord * 0.1f);
                streamline9.transform.position = new UnityEngine.Vector3(0.25f, 0.25f, 0);
                streamline9.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                streamline9.mesh.vertices = newVertices;
                streamline9.mesh.SetIndices(indicesForStreamlines, MeshTopology.LineStrip, 0, true);
                streamline9.GetComponent<Renderer>().enabled = true;
            }
        }
        else if (solver.airfoilType != solver.DisplayGrid.joukowsky || solver.airfoilType == solver.DisplayGrid.joukowsky && solver.displayJoukowsky == solver.DisplayJoukowsky.numericSolver)
        {
            airfoilType = -1;
            streamline0.GetComponent<Renderer>().enabled = false;
            streamline1.GetComponent<Renderer>().enabled = false;
            streamline2.GetComponent<Renderer>().enabled = false;
            streamline3.GetComponent<Renderer>().enabled = false;
            streamline4.GetComponent<Renderer>().enabled = false;
            streamline5.GetComponent<Renderer>().enabled = false;
            streamline6.GetComponent<Renderer>().enabled = false;
            streamline7.GetComponent<Renderer>().enabled = false;
            streamline8.GetComponent<Renderer>().enabled = false;
            streamline9.GetComponent<Renderer>().enabled = false;
        }
    }

    void Start()
    {
        vertices = streamline0.mesh.vertices = streamline1.mesh.vertices = streamline2.mesh.vertices = streamline3.mesh.vertices = streamline4.mesh.vertices = streamline5.mesh.vertices = streamline6.mesh.vertices = streamline7.mesh.vertices = streamline8.mesh.vertices = streamline9.mesh.vertices;
        newVertices = new UnityEngine.Vector3[101];
        zeroVertices = new UnityEngine.Vector3[101];
        for (int i = 0; i < 101; ++i)
        {
            zeroVertices[i] = UnityEngine.Vector3.zero;
            zeroIndices[i] = 0;
        }
    }

    void Update()
    {
        Joukowsky();
    }
}
