using System;
using System.Numerics;
using UnityEngine;

public class solver : MonoBehaviour
{
    MeshFilter airfoil;
    BoxCollider boxCollider;
    [HideInInspector] public float chord;
    [Range(0, 20f)] public float angleOfAttack;
    [HideInInspector] public float angleOfAttack0 = 0;
    [Space]
    [Header("grid")]
    public int resolution = 100;
    public DisplayGrid airfoilType;
    [HideInInspector] public DisplayGrid airfoilType0 = 0;
    public enum DisplayGrid { joukowsky, nacaFourDigit, nacaFourDigitModified, nacaFiveDigit, nacaFiveDigitModified };
    [Space]
    [Header("inlet emission")]
    public float injectInletDensity = 1;
    public float inletRadius = 0;
    [Range(0, 100)] public int inletSamples;
    [Space]
    [Header("mouse emission")]
    public float injectDensity = 1;
    public float injectVelocity = 1;
    public float maxInjectVelocity = 10;
    public float radius = 0;
    [Range(1, 100)] public int samples = 1;
    [Space]
    [Header("simulation")]
    [Range(0, 1)] public float friction = 0.001f;
    [Range(0, 10)] public float acceleration = 1;
    [Range(0, 1)] public float swirl = 0;
    [Space]
    [Header("diffusion")]
    public bool diffuseDensity = true;
    public bool diffuseVelocity = true;
    [Range(0, 1)] public float densityDiffusion = 1;
    [Range(0, 1)] public float velocityDiffusion = 0.1f;
    [Range(0, 5)] public int diffusionQuality = 0;
    [Space]
    [Header("advection")]
    public bool advectDensity = true;
    public enum VelocityAdvection { constant, variable };
    public VelocityAdvection velocityAdvectionType;
    public float densityAdvection = 1;
    public float velocityAdvection = 1;
    [Space]
    [Header("dissipation")]
    public float densityDissipation = 0;
    public float velocityDissipation = 0;
    [Space]
    [Header("conservation")]
    public bool conserveMass = true;
    [Range(1, 50)] public int solverQuality = 20;
    [Space]
    [Header("visualization")]
    public bool displayVelocity;
    public float displayVelocityMultiplier = 1;
    public enum DisplayG { none, density, pressure };
    public DisplayG displayG;
    public float displayGridMultiplier = 1;
    [Space]
    [Header("joukowsky")]
    public DisplayJoukowsky displayJoukowsky;
    [HideInInspector] public DisplayJoukowsky displayJoukowsky0 = 0;
    public enum DisplayJoukowsky { analyticVelocity, analyticPressure, numericSolver };
    [Range(0, 0.25f)] public float xDisplacement;
    [HideInInspector] public float xDisplacement0 = 0;
    [Range(0, 0.25f)] public float yDisplacement;
    [HideInInspector] public float yDisplacement0 = 0;
    public float freestreamVelocity;
    [HideInInspector] public float freestreamVelocity0 = 0;
    [Space]
    [Header("naca 4-digit airfoil")]
    public float maxCamber4;
    float maxCamber40;
    public float maxCamberPosition4;
    float maxCamberPosition40 = 0;
    public float maxThickness4 = 0;
    float maxThickness40 = 0;
    [Space]
    [Header("naca 4-digit modified airfoil")]
    public float maxCamber4m;
    float maxCamber4m0 = 0;
    public float maxCamberPosition4m;
    float maxCamberPosition4m0 = 0;
    public float maxThickness4m;
    float maxThickness4m0 = 0;
    public enum Display4m { SixtyTwo = 62, SixtyThree = 63, SixtyFour = 64, SixtyFive = 65, Three = 3, ThirtyThree = 33, NinetyThree = 93, Five = 5, ThirtyFive = 35, ThirtyFour = 34 };
    public Display4m modifications4m;
    [HideInInspector] public Display4m modifications4m0 = 0;
    [Space]
    [Header("naca 5-digit airfoil")]
    public float designLiftCoefficient5;
    float designLiftCoefficient50 = 0;
    [Range(1, 5)] public int maxCamberPosition5 = 1;
    int maxCamberPosition50 = 1;
    [Range(0, 1)] public int camberType5 = 0;
    int camberType50 = 0;
    public float maxThickness5;
    float maxThickness50 = 0;
    [Space]
    [Header("naca 5-digit modified airfoil")]
    public float designLiftCoefficient5m;
    float designLiftCoefficient5m0 = 0;
    [Range(1, 5)] public int maxCamberPosition5m = 1;
    int maxCamberPosition5m0 = 1;
    [Range(0, 1)] public int camberType5m = 0;
    int camberType5m0 = 0;
    public float maxThickness5m;
    float maxThickness5m0 = 0;
    public enum Display5m { SixtyTwo = 62, SixtyThree = 63, SixtyFour = 64, SixtyFive = 65, Three = 3, ThirtyThree = 33, NinetyThree = 93, Five = 5, ThirtyFive = 35, ThirtyFour = 34 };
    public Display5m modifications5m;
    [HideInInspector] public Display5m modifications5m0 = 0;
    float maxCamberActual;
    float maxCamberPositionActual;
    float thicknessActual;
    float designLiftCoefficientActual;
    UnityEngine.Vector3[] vertices, circle;
    float chordActual;
    float position;
    float camberFunction;
    float derivativeCamberFunction;
    float angle;
    float thickness;
    float r, k1, k1k2, a0, a1, a2, a3, d0, d1, d2, d3;

    int[] indicesForStreamlines;
    int voxelCount;
    float voxelScale, voxelScale2;
    int resolution0, resolution2 = 0;
    float[] velocityMagnitude, joukowskyPressure, density0, density, divergence, pressure, vorticity, absoluteVorticity, vorticityLength;
    UnityEngine.Vector3 v, v1, scale;
    UnityEngine.Vector2[] velocity0, velocity, vorticityGradient, vorticityConfinement;
    Color[] colors = new Color[101];

    ParticleSystem.Particle[] particles_scalarField, particles_vectorField;
    ParticleSystem particleSystem_scalarField, particleSystem_vectorField;
    Renderer renderer_scalarField, renderer_vectorField;
    ParticleCollisionEvent[] collisionEvents;

    void Initialize()
    {
        if (resolution0 != resolution)
        {
            scale = 2f * transform.localScale;
            if (scale.x == 0.0f || resolution <= 0) return;

            resolution0 = resolution;
            resolution2 = resolution + 2;
            voxelScale = scale.x / resolution;
            voxelScale2 = scale.x / resolution2;
            voxelCount = resolution2 * resolution2;

            velocityMagnitude = joukowskyPressure = new float[voxelCount];
            indicesForStreamlines = new int[101];
            density0 = density = new float[voxelCount];
            velocity0 = velocity = new UnityEngine.Vector2[voxelCount];
            divergence = new float[voxelCount];
            pressure = new float[voxelCount];
            vorticity = new float[voxelCount];
            absoluteVorticity = new float[voxelCount];
            vorticityLength = new float[voxelCount];
            vorticityGradient = new UnityEngine.Vector2[voxelCount];
            vorticityConfinement = new UnityEngine.Vector2[voxelCount];

            int ii = 0;
            particles_scalarField = new ParticleSystem.Particle[resolution * resolution];
            particles_vectorField = new ParticleSystem.Particle[resolution * resolution];

            for (int i = resolution / 3; i < (2 * resolution) / 3; ++i)
                for (int j = 0; j < resolution; ++j)
                {
                    v.x = i * voxelScale;
                    v.y = j * voxelScale;
                    particles_scalarField[ii].position = v;
                    particles_scalarField[ii].startSize = voxelScale;
                    particles_vectorField[ii].position = v + 0.9f * UnityEngine.Vector3.left;
                    particles_vectorField[ii].startSize = 0.001f;
                    collisionEvents = new ParticleCollisionEvent[ii];
                    ++ii;
                }
            airfoilType0 = airfoilType;
            displayJoukowsky0 = displayJoukowsky;
            var emitParams = new ParticleSystem.EmitParams();
            emitParams.position = UnityEngine.Vector3.zero;
            particleSystem_scalarField.Emit(emitParams, resolution * resolution);
            particleSystem_scalarField.SetParticles(particles_scalarField, resolution * resolution);
            particleSystem_vectorField.Emit(emitParams, resolution * resolution);
            particleSystem_vectorField.SetParticles(particles_vectorField, resolution * resolution);
        }
    }

    UnityEngine.Vector2 position0 = new UnityEngine.Vector2();
    UnityEngine.Vector2 bend = new UnityEngine.Vector2();

    //experiment with velocity emission hit points, position, v, etc.
    void Emit()
    {
        if (Input.GetMouseButton(0) || Input.GetMouseButton(1))
        {
            if (Physics.Raycast(Camera.main.ScreenPointToRay(Input.mousePosition), out RaycastHit hit))
            {
                // Emit multiple times (defined by the number of "samples") within the cursor area (defined by its "radius" parameter).
                int idx;
                float s = 2f * transform.localScale.x;
                for (int i = 0; i < samples; ++i)
                {
                    if (samples == 1) v = UnityEngine.Vector2.zero;
                    else v = UnityEngine.Quaternion.AngleAxis(UnityEngine.Random.Range(0, 360), UnityEngine.Vector3.forward) * UnityEngine.Vector3.right * UnityEngine.Random.Range(0, radius);
                    v += hit.point;
                    // Skip if the emission point is outside the voxel grid.
                    if (v.x < 0 || v.x > s) return;
                    if (v.y < 0 || v.y > s) return;
                    // Calculate the voxel index.
                    v.x = (float)Math.Ceiling(v.x / voxelScale);
                    v.y = (float)Math.Ceiling((v.y - 0.66f) / voxelScale);
                    idx = (int)(v.y * resolution2 + v.x);
                    if (idx < 0) return;

                    // Add density.
                    if (Input.GetMouseButton(0))//(Input.GetKey(KeyCode.Space))
                    {
                        density0[idx] = injectDensity;
                    }

                    // Add velocity.
                    if (Input.GetMouseButtonDown(1))
                    {
                        position0.x = hit.point.x;
                        position0.y = hit.point.y;
                    }
                    else if (Input.GetMouseButton(1))
                    {
                        // Calculate direction vector from previous and current cursor positions.
                        v.x = position0.x = hit.point.x - position0.x;
                        v.y = position0.y = hit.point.y - position0.y;
                        // Bend it towards the mouse direction, so the result feels more natural to the user.
                        position0.x += bend.x;
                        position0.y += bend.y;
                        velocity0[idx] += position0.normalized * injectVelocity;
                        // Clamp the result, so it does not overshoot.
                        if (velocity0[idx].magnitude > maxInjectVelocity)
                            velocity0[idx] = velocity0[idx].normalized * maxInjectVelocity;
                        // Remember the current bend vector.
                        bend.x = v.x;
                        bend.y = v.y;
                        // Remember the current cursor position.
                        position0.x = hit.point.x;
                        position0.y = hit.point.y;
                    }
                }
            }
        }

        if (inletSamples != 0)
        {
            int idx0;
            float s = 2f * transform.localScale.x;
            for (int i = 0; i < inletSamples; ++i)
            {
                if (inletSamples == 1) v1 = UnityEngine.Vector2.zero;
                else v1 = UnityEngine.Quaternion.AngleAxis(UnityEngine.Random.Range(0, 360), UnityEngine.Vector3.forward) * UnityEngine.Vector3.right * UnityEngine.Random.Range(0, inletRadius);
                v1 += new UnityEngine.Vector3(0.25f, 1f);
                // Calculate the voxel index.
                v1.x = (float)Math.Ceiling(v1.x / voxelScale);
                v1.y = (float)Math.Ceiling((v1.y - 0.66f) / voxelScale);
                idx0 = (int)(v1.y * resolution2 + v1.x);
                if (idx0 < 0 || idx0 > voxelCount) return;

                // Add density.
                density0[idx0] = injectInletDensity;
            }
        }
    }

    // Generate small scale vorticles to increase detail in the simulation.
    void VorticityConfinement()
    {
        if (swirl <= 0) return; else if (swirl > 1) swirl = 1;

        int ii;
        for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                ii = i + resolution2 * j;
                vorticity[ii] = (velocity0[ii + resolution2].y - velocity0[ii - resolution2].y - velocity0[ii + 1].x + velocity0[ii - 1].x) * 0.5f;
                if (vorticity[ii] >= 0) absoluteVorticity[ii] = vorticity[ii];
                else absoluteVorticity[ii] = -vorticity[ii];
            }

        UnityEngine.Vector2 v;
        for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                ii = i + resolution2 * j;
                v.x = (absoluteVorticity[ii + 1] - absoluteVorticity[ii - 1]) * 0.5f;
                v.y = (absoluteVorticity[ii + resolution2] - absoluteVorticity[ii - resolution2]) * 0.5f;
                vorticityGradient[ii] = v;
                vorticityLength[ii] = Mathf.Sqrt(vorticityGradient[ii].x * vorticityGradient[ii].x + vorticityGradient[ii].y * vorticityGradient[ii].y);
                if (vorticityLength[ii] < 0.01f)
                    vorticityConfinement[ii] = UnityEngine.Vector2.zero;
                else
                {
                    v.x = vorticityGradient[ii].x / vorticityLength[ii];
                    v.x = vorticityGradient[ii].y / vorticityLength[ii];
                    vorticityConfinement[ii] = v;
                }
            }

        for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                ii = i + resolution2 * j;
                v.x = swirl * (vorticityConfinement[ii].y * vorticity[ii]);
                v.y = swirl * (-vorticityConfinement[ii].x * vorticity[ii]);
                velocity0[ii] += v;
            }
    }

    // Diffuse the density/velocity fields.
    void Diffuse()
    {
        // Make sure that dissipation is computed (if the related parameters are non-zero) regardless of the state of the diffuse density/velocity flags.
        if (!diffuseDensity && densityDissipation > 0) diffuseDensity = true;
        if (!diffuseVelocity && velocityDissipation > 0) diffuseVelocity = true;

        if (!diffuseDensity && !diffuseVelocity) return;

        if (densityDissipation < 0) densityDissipation = 0;
        if (velocityDissipation < 0) velocityDissipation = 0;

        int ii;
        float idd = 1 - densityDiffusion;
        float idv = 1 - velocityDiffusion;
        float a = densityDiffusion * resolution * resolution * Time.deltaTime;
        float b = velocityDiffusion * resolution * resolution * Time.deltaTime;
        for (int jj = 0; jj <= diffusionQuality; ++jj)
            for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    ii = i + resolution2 * j;
                    if (diffuseDensity)
                    {
                        density[ii] = (density0[ii] * idd + (((density0[ii] + a * (density0[ii - 1] + density0[ii + 1] + density0[ii - resolution2] + density0[ii + resolution2])) / (1 + 4 * a)) * densityDiffusion)) - densityDissipation;
                        if (density[ii] < 0) density[ii] = 0;
                    }
                    if (diffuseVelocity)
                    {
                        velocity[ii] = (velocity0[ii] * idv + (((velocity0[ii] + b * (velocity0[ii - 1] + velocity0[ii + 1] + velocity0[ii - resolution2] + velocity0[ii + resolution2])) / (1 + 4 * b)) * velocityDiffusion));
                        if (velocityDissipation > 0)
                        {
                            float m = velocity[ii].magnitude - velocityDissipation;
                            if (m < 0) m = 0;
                            velocity[ii] = velocity[ii].normalized * m;
                        }
                    }
                }
        if (diffuseDensity) density0 = density;
        if (diffuseVelocity) velocity0 = velocity;
    }

    // Advecty densitiy along the velocity field using backward advection.
    void AdvectDensity()
    {
        if (!advectDensity) return;

        int ii, i0, j0, i1, j1;
        float x, y, s1, s0, t1, t0;
        float c = Time.deltaTime * resolution * densityAdvection;
        density = new float[voxelCount];
        for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                ii = i + resolution2 * j;
                x = i - velocity0[ii].x * c;
                y = j - velocity0[ii].y * c;
                if (x < 0.5f) x = 0.5f; else if (x > resolution + 0.5f) x = resolution + 0.5f; i0 = (int)x; i1 = i0 + 1;
                if (y < 0.5f) y = 0.5f; else if (y > resolution + 0.5f) y = resolution + 0.5f; j0 = (int)y; j1 = j0 + 1;
                s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
                density[ii] = s0 * (t0 * density0[i0 + resolution2 * j0] + t1 * density0[i0 + resolution2 * j1]) + s1 * (t0 * density0[i1 + resolution2 * j0] + t1 * density0[i1 + resolution2 * j1]);
            }
        density0 = density;
    }

    // Self-advect velocity field using backward advection.
    void AdvectVelocity()
    {
        if (velocityAdvectionType == VelocityAdvection.variable)
        {
            if (friction <= 0) friction = 0; else if (friction > 1) friction = 1;

            int ii, i0, j0, i1, j1;
            float x, y, s1, s0, t1, t0;
            float c = Time.deltaTime * resolution * velocityAdvection;
            float f = 1 - friction;
            UnityEngine.Vector2 g = new UnityEngine.Vector2(acceleration * 0.00001f, 0);
            velocity = new UnityEngine.Vector2[voxelCount];
            for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    ii = i + resolution2 * j;
                    x = i - velocity0[ii].x * c;
                    y = j - velocity0[ii].y * c;
                    if (x < 0.5f) x = 0.5f; else if (x > resolution + 0.5f) x = resolution + 0.5f; i0 = (int)x; i1 = i0 + 1;
                    if (y < 0.5f) y = 0.5f; else if (y > resolution + 0.5f) y = resolution + 0.5f; j0 = (int)y; j1 = j0 + 1;
                    s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
                    velocity[ii] = (s0 * (t0 * velocity0[i0 + resolution2 * j0] + t1 * velocity0[i0 + resolution2 * j1]) + s1 * (t0 * velocity0[i1 + resolution2 * j0] + t1 * velocity0[i1 + resolution2 * j1])) * f + density0[ii] * g;
                }
            velocity0 = velocity;
        }

        else if (velocityAdvectionType == VelocityAdvection.constant)
        {
            int ii;
            for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    ii = i + resolution2 * j;
                    velocity[ii] = 0.25f * UnityEngine.Vector2.right;
                }
            velocity0 = velocity;
        }
    }

    // Mass conservation.
    void Project()
    {
        if (!conserveMass) return;

        int ii;
        float h = 1.0f / resolution;
        for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                ii = i + resolution2 * j;
                divergence[ii] = -0.5f * h * (velocity0[ii + 1].x - velocity0[ii - 1].x + velocity0[ii + resolution2].y - velocity0[ii - resolution2].y);
            }

        for (int k = 0; k < solverQuality; ++k)
            for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    ii = i + resolution2 * j;
                    pressure[ii] = (divergence[ii] + pressure[ii - 1] + pressure[ii + 1] + pressure[ii - resolution2] + pressure[ii + resolution2]) / 4;
                }

        for (int i = 1; i <= resolution; ++i)
            for (int j = 1; j <= resolution; ++j)
            {
                ii = i + resolution2 * j;
                velocity0[ii].x -= 0.5f * (pressure[ii + 1] - pressure[ii - 1]) / h;
                velocity0[ii].y -= 0.5f * (pressure[ii + resolution2] - pressure[ii - resolution2]) / h;
            }
    }

    void Display()
    {
        // Scalar fields.
        renderer_scalarField.enabled = true;
        if (airfoilType == DisplayGrid.joukowsky)
        {
            float xDisplacement00 = 0, yDisplacement00 = 0, angleOfAttack00 = 0, freestreamVelocity00 = 0;
            if (xDisplacement00 != xDisplacement || yDisplacement00 != yDisplacement || angleOfAttack00 != angleOfAttack || freestreamVelocity00 != freestreamVelocity)
            {
                xDisplacement00 = xDisplacement;
                yDisplacement00 = yDisplacement;
                angleOfAttack00 = angleOfAttack;
                freestreamVelocity00 = freestreamVelocity;
                int ii = 0;
                float value;

                if (displayJoukowsky == DisplayJoukowsky.analyticVelocity)
                {
                    for (int i = 1; i <= resolution; ++i)
                    {
                        for (int j = 1; j <= resolution; ++j)
                        {
                            value = 0.02f * velocityMagnitude[j + resolution2 * i];
                            particles_scalarField[ii].startColor = new Color(value, value, value, 1f);
                            ++ii;
                        }
                    }
                }

                else if (displayJoukowsky == DisplayJoukowsky.analyticPressure)
                {
                    for (int i = 1; i <= resolution; ++i)
                    {
                        for (int j = 1; j <= resolution; ++j)
                        {
                            value = 0.15f * (4f + joukowskyPressure[j + resolution2 * i]);
                            particles_scalarField[ii].startColor = new Color(value, value, value, 1f);
                            ++ii;
                        }
                    }
                }

                else if (displayJoukowsky == DisplayJoukowsky.numericSolver)
                {
                    int jj, min_X, max_X, min_Y, max_Y, min_XX, min_XY, max_XX, max_XY, min_YX, min_YY, max_YX, max_YY;
                    int[] minBoundaryIndices, maxBoundaryIndices;
                    float alpha;
                    UnityEngine.Vector3 minX, maxX, minY, maxY;
                    int[] boundaryIndices = new int[vxlcnt], boundaryIndices_X = new int[vxlcnt], boundaryIndices_Y = new int[vxlcnt];
                    renderer_scalarField.enabled = true;

                    if (displayG == DisplayG.density)
                    {
                        ii = jj = 0;
                        minX = maxX = minY = maxY = vertices[0];
                        alpha = angleOfAttack * Mathf.PI / 180f;

                        for (int i = 0; i < vxlcnt; ++i)
                        {
                            if (minX.x > vertices[i].x) minX = vertices[i];
                            if (maxX.x < vertices[i].x) maxX = vertices[i];
                            if (minY.y > vertices[i].y) minY = vertices[i];
                            if (maxY.y < vertices[i].y) maxY = vertices[i];
                            boundaryIndices[i] = PosToIdx(transform.TransformPoint(0.5f * (vertices[i] - new UnityEngine.Vector3(0.75f, 1f))));
                            boundaryIndices_Y[i] = (int)Math.Ceiling((double)(boundaryIndices[i] / resolution));
                            boundaryIndices_X[i] = boundaryIndices[i] - boundaryIndices_Y[i] * resolution;
                        }

                        min_X = PosToIdx(transform.TransformPoint(0.5f * (minX - new UnityEngine.Vector3(0.75f, 1f))));
                        max_X = PosToIdx(transform.TransformPoint(0.5f * (maxX - new UnityEngine.Vector3(0.75f, 1f))));
                        min_Y = PosToIdx(transform.TransformPoint(0.5f * (minY - new UnityEngine.Vector3(0.75f, 1f))));
                        max_Y = PosToIdx(transform.TransformPoint(0.5f * (maxY - new UnityEngine.Vector3(0.75f, 1f))));

                        min_XY = (int)Math.Ceiling((double)(min_X / resolution));
                        min_XX = min_X - min_XY * resolution;
                        max_XY = (int)Math.Ceiling((double)(max_X / resolution));
                        max_XX = max_X - max_XY * resolution;
                        min_YY = (int)Math.Ceiling((double)(min_Y / resolution));
                        min_YX = min_Y - min_YY * resolution;
                        max_YY = (int)Math.Ceiling((double)(max_Y / resolution));
                        max_YX = max_Y - max_YY * resolution;

                        minBoundaryIndices = new int[max_XX - min_XX + 1];
                        maxBoundaryIndices = new int[max_XX - min_XX + 1];

                        for (int i = min_XX; i <= max_XX; ++i)
                        {
                            minBoundaryIndices[i - min_XX] = resolution2 * resolution2 + 1;
                            maxBoundaryIndices[i - min_XX] = -1;
                        }

                        for (int i = 1; i <= resolution; ++i)
                        {
                            for (int j = 1; j <= resolution; ++j)
                            {
                                value = density[j + resolution2 * i] * displayGridMultiplier;
                                particles_scalarField[ii].startColor = new Color(value, value, value, 1f);
                                ++ii;
                            }
                        }

                        for (int i = 0; i < vxlcnt; ++i)
                        {
                            for (int k = min_XX; k <= max_XX; ++k)
                            {
                                for (int j = min_YY; j <= max_YY; ++j)
                                {
                                    jj = k + resolution * j;
                                    if (boundaryIndices_X[i] == k)
                                    {
                                        if (minBoundaryIndices[k - min_XX] > boundaryIndices[i]) minBoundaryIndices[k - min_XX] = boundaryIndices[i];
                                        if (maxBoundaryIndices[k - min_XX] < boundaryIndices[i]) maxBoundaryIndices[k - min_XX] = boundaryIndices[i];
                                        if (jj >= minBoundaryIndices[k - min_XX] && jj <= maxBoundaryIndices[k - min_XX])
                                        {
                                            density0[k - (int)Mathf.Floor(resolution * 0.125f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = 0;
                                            velocity0[k - (int)Mathf.Floor(resolution * 0.125f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = UnityEngine.Vector2.zero;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    else if (displayG == DisplayG.pressure)
                    {
                        velocityAdvectionType = VelocityAdvection.variable;
                        ii = jj = 0;

                        minX = maxX = minY = maxY = vertices[0];
                        alpha = angleOfAttack * Mathf.PI / 180f;

                        for (int i = 0; i < vxlcnt; ++i)
                        {
                            if (minX.x > vertices[i].x) minX = vertices[i];
                            if (maxX.x < vertices[i].x) maxX = vertices[i];
                            if (minY.y > vertices[i].y) minY = vertices[i];
                            if (maxY.y < vertices[i].y) maxY = vertices[i];
                            boundaryIndices[i] = PosToIdx(transform.TransformPoint(0.5f * (vertices[i] - new UnityEngine.Vector3(0.75f, 1f))));
                            boundaryIndices_Y[i] = (int)Math.Ceiling((double)(boundaryIndices[i] / resolution));
                            boundaryIndices_X[i] = boundaryIndices[i] - boundaryIndices_Y[i] * resolution;
                        }

                        min_X = PosToIdx(transform.TransformPoint(0.5f * (minX - new UnityEngine.Vector3(0.75f, 1f))));
                        max_X = PosToIdx(transform.TransformPoint(0.5f * (maxX - new UnityEngine.Vector3(0.75f, 1f))));
                        min_Y = PosToIdx(transform.TransformPoint(0.5f * (minY - new UnityEngine.Vector3(0.75f, 1f))));
                        max_Y = PosToIdx(transform.TransformPoint(0.5f * (maxY - new UnityEngine.Vector3(0.75f, 1f))));

                        min_XY = (int)Math.Ceiling((double)(min_X / resolution));
                        min_XX = min_X - min_XY * resolution;
                        max_XY = (int)Math.Ceiling((double)(max_X / resolution));
                        max_XX = max_X - max_XY * resolution;
                        min_YY = (int)Math.Ceiling((double)(min_Y / resolution));
                        min_YX = min_Y - min_YY * resolution;
                        max_YY = (int)Math.Ceiling((double)(max_Y / resolution));
                        max_YX = max_Y - max_YY * resolution;

                        minBoundaryIndices = new int[max_YY - min_YY + 1];
                        maxBoundaryIndices = new int[max_YY - min_YY + 1];

                        for (int i = min_YY; i <= max_YY; ++i)
                        {
                            minBoundaryIndices[i - min_YY] = resolution2 * resolution2 + 1;
                            maxBoundaryIndices[i - min_YY] = -1;
                        }

                        for (int i = 0; i < vxlcnt; ++i)
                        {
                            for (int j = min_YY; j <= max_YY; ++j)
                            {
                                for (int k = min_XX; k <= max_XX; ++k)
                                {
                                    jj = k + resolution * j;
                                    if (boundaryIndices_Y[i] == j)
                                    {
                                        if (minBoundaryIndices[j - min_YY] > boundaryIndices[i]) minBoundaryIndices[j - min_YY] = boundaryIndices[i];
                                        if (maxBoundaryIndices[j - min_YY] < boundaryIndices[i]) maxBoundaryIndices[j - min_YY] = boundaryIndices[i];
                                        if (jj > minBoundaryIndices[j - min_YY] && jj < maxBoundaryIndices[j - min_YY])
                                        {
                                            pressure[k - (int)Mathf.Floor(resolution * 0.125f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = 0;
                                        }
                                        else if (jj == minBoundaryIndices[j - min_YY] || jj == maxBoundaryIndices[j - min_YY])
                                        {
                                            pressure[k - (int)Mathf.Floor(resolution * 0.125f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = 0.05f * Mathf.Cos(2f * alpha) - Mathf.Pow(velocity0[jj].x, 2);
                                        }
                                    }
                                }
                            }
                        }

                        for (int i = 1; i <= resolution; ++i)
                            for (int j = 1; j <= resolution; ++j)
                            {
                                value = pressure[j + resolution2 * i] * displayGridMultiplier * 50000;
                                particles_scalarField[ii].startColor = new Color(value, value, value, 1);
                                ++ii;
                            }
                    }
                }
            }
            particleSystem_scalarField.SetParticles(particles_scalarField, voxelCount);
        }

        else if (airfoilType != DisplayGrid.joukowsky && displayG != DisplayG.none)
        {
            int ii, jj, min_X, max_X, min_Y, max_Y, min_XX, min_XY, max_XX, max_XY, min_YX, min_YY, max_YX, max_YY;
            int[] minBoundaryIndices, maxBoundaryIndices;
            float value, alpha;
            UnityEngine.Vector3 minX, maxX, minY, maxY;
            int[] boundaryIndices = new int[vxlcnt], boundaryIndices_X = new int[vxlcnt], boundaryIndices_Y = new int[vxlcnt];
            renderer_scalarField.enabled = true;
            if (displayG == DisplayG.density)
            {
                ii = jj = 0;
                minX = maxX = minY = maxY = vertices[0];
                alpha = angleOfAttack * Mathf.PI / 180f;

                for (int i = 0; i < vxlcnt; ++i)
                {
                    if (minX.x > vertices[i].x) minX = vertices[i];
                    if (maxX.x < vertices[i].x) maxX = vertices[i];
                    if (minY.y > vertices[i].y) minY = vertices[i];
                    if (maxY.y < vertices[i].y) maxY = vertices[i];
                    boundaryIndices[i] = PosToIdx(transform.TransformPoint(0.5f * (vertices[i] - new UnityEngine.Vector3(0.75f, 1f))));
                    boundaryIndices_Y[i] = (int)Math.Ceiling((double)(boundaryIndices[i] / resolution));
                    boundaryIndices_X[i] = boundaryIndices[i] - boundaryIndices_Y[i] * resolution;
                }

                min_X = PosToIdx(transform.TransformPoint(0.5f * (minX - new UnityEngine.Vector3(0.75f, 1f))));
                max_X = PosToIdx(transform.TransformPoint(0.5f * (maxX - new UnityEngine.Vector3(0.75f, 1f))));
                min_Y = PosToIdx(transform.TransformPoint(0.5f * (minY - new UnityEngine.Vector3(0.75f, 1f))));
                max_Y = PosToIdx(transform.TransformPoint(0.5f * (maxY - new UnityEngine.Vector3(0.75f, 1f))));

                min_XY = (int)Math.Ceiling((double)(min_X / resolution));
                min_XX = min_X - min_XY * resolution;
                max_XY = (int)Math.Ceiling((double)(max_X / resolution));
                max_XX = max_X - max_XY * resolution;
                min_YY = (int)Math.Ceiling((double)(min_Y / resolution));
                min_YX = min_Y - min_YY * resolution;
                max_YY = (int)Math.Ceiling((double)(max_Y / resolution));
                max_YX = max_Y - max_YY * resolution;

                minBoundaryIndices = new int[max_XX - min_XX + 1];
                maxBoundaryIndices = new int[max_XX - min_XX + 1];

                for (int i = min_XX; i <= max_XX; ++i)
                {
                    minBoundaryIndices[i - min_XX] = resolution2 * resolution2 + 1;
                    maxBoundaryIndices[i - min_XX] = -1;
                }

                for (int i = 1; i <= resolution; ++i)
                {
                    for (int j = 1; j <= resolution; ++j)
                    {
                        value = density[j + resolution2 * i] * displayGridMultiplier;
                        particles_scalarField[ii].startColor = new Color(value, value, value, 1f);
                        ++ii;
                    }
                }

                for (int i = 0; i < vxlcnt; ++i)
                {
                    for (int k = min_XX; k <= max_XX; ++k)
                    {
                        for (int j = min_YY; j <= max_YY; ++j)
                        {
                            jj = k + resolution * j;
                            if (boundaryIndices_X[i] == k)
                            {
                                if (minBoundaryIndices[k - min_XX] > boundaryIndices[i]) minBoundaryIndices[k - min_XX] = boundaryIndices[i];
                                if (maxBoundaryIndices[k - min_XX] < boundaryIndices[i]) maxBoundaryIndices[k - min_XX] = boundaryIndices[i];
                                if (jj >= minBoundaryIndices[k - min_XX] && jj <= maxBoundaryIndices[k - min_XX])
                                {
                                    density0[k + 2 - (int)Mathf.Floor(resolution * 0.01f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = 0;
                                    velocity0[k + 2 - (int)Mathf.Floor(resolution * 0.01f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = UnityEngine.Vector2.zero;
                                }
                            }
                        }
                    }
                }
            }

            else if (displayG == DisplayG.pressure)
            {
                velocityAdvectionType = VelocityAdvection.variable;
                ii = jj = 0;

                minX = maxX = minY = maxY = vertices[0];
                alpha = angleOfAttack * Mathf.PI / 180f;

                for (int i = 0; i < vxlcnt; ++i)
                {
                    if (minX.x > vertices[i].x) minX = vertices[i];
                    if (maxX.x < vertices[i].x) maxX = vertices[i];
                    if (minY.y > vertices[i].y) minY = vertices[i];
                    if (maxY.y < vertices[i].y) maxY = vertices[i];
                    boundaryIndices[i] = PosToIdx(transform.TransformPoint(0.5f * (vertices[i] - new UnityEngine.Vector3(0.75f, 1f))));
                    boundaryIndices_Y[i] = (int)Math.Ceiling((double)(boundaryIndices[i] / resolution));
                    boundaryIndices_X[i] = boundaryIndices[i] - boundaryIndices_Y[i] * resolution;
                }

                min_X = PosToIdx(transform.TransformPoint(0.5f * (minX - new UnityEngine.Vector3(0.75f, 1f))));
                max_X = PosToIdx(transform.TransformPoint(0.5f * (maxX - new UnityEngine.Vector3(0.75f, 1f))));
                min_Y = PosToIdx(transform.TransformPoint(0.5f * (minY - new UnityEngine.Vector3(0.75f, 1f))));
                max_Y = PosToIdx(transform.TransformPoint(0.5f * (maxY - new UnityEngine.Vector3(0.75f, 1f))));

                min_XY = (int)Math.Ceiling((double)(min_X / resolution));
                min_XX = min_X - min_XY * resolution;
                max_XY = (int)Math.Ceiling((double)(max_X / resolution));
                max_XX = max_X - max_XY * resolution;
                min_YY = (int)Math.Ceiling((double)(min_Y / resolution));
                min_YX = min_Y - min_YY * resolution;
                max_YY = (int)Math.Ceiling((double)(max_Y / resolution));
                max_YX = max_Y - max_YY * resolution;

                minBoundaryIndices = new int[max_YY - min_YY + 1];
                maxBoundaryIndices = new int[max_YY - min_YY + 1];

                for (int i = min_YY; i <= max_YY; ++i)
                {
                    minBoundaryIndices[i - min_YY] = resolution2 * resolution2 + 1;
                    maxBoundaryIndices[i - min_YY] = -1;
                }

                for (int i = 0; i < vxlcnt; ++i)
                {
                    for (int k = min_XX; k <= max_XX; ++k)
                    {
                        for (int j = min_YY; j <= max_YY; ++j)
                        {
                            jj = k + resolution * j;
                            if (boundaryIndices[i] > j * resolution && boundaryIndices[i] < (j + 1) * resolution)
                            {
                                if (minBoundaryIndices[j - min_YY] > boundaryIndices[i]) minBoundaryIndices[j - min_YY] = boundaryIndices[i];
                                if (maxBoundaryIndices[j - min_YY] < boundaryIndices[i]) maxBoundaryIndices[j - min_YY] = boundaryIndices[i];
                                if (jj > minBoundaryIndices[j - min_YY] && jj < maxBoundaryIndices[j - min_YY])
                                {
                                    pressure[k + 2 - (int)Mathf.Floor(resolution * 0.01f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = 0;
                                }
                                if (jj == minBoundaryIndices[j - min_YY] || jj == maxBoundaryIndices[j - min_YY])
                                {
                                    pressure[k + 2 - (int)Mathf.Floor(resolution * 0.01f) + resolution2 * (j - (int)Mathf.Floor(resolution * 0.325f))] = 0.05f * Mathf.Cos(2f * alpha) - Mathf.Pow(velocity0[jj].x, 2);
                                }
                            }
                        }
                    }
                }

                for (int i = 1; i <= resolution; ++i)
                    for (int j = 1; j <= resolution; ++j)
                    {
                        value = pressure[j + resolution2 * i] * displayGridMultiplier * 50000;
                        particles_scalarField[ii].startColor = new Color(value, value, value, 1f);
                        ++ii;
                    }
            }
            particleSystem_scalarField.SetParticles(particles_scalarField, voxelCount);
        }
        else renderer_scalarField.enabled = false;

        if (displayVelocity && airfoilType != DisplayGrid.joukowsky || displayVelocity && airfoilType == DisplayGrid.joukowsky && displayJoukowsky == DisplayJoukowsky.numericSolver)
        {
            int ii = 0;
            renderer_vectorField.enabled = true;
            for (int i = 1; i <= resolution; ++i)
                for (int j = 1; j <= resolution; ++j)
                {
                    particles_vectorField[ii].velocity = -new UnityEngine.Vector2(velocity0[j + resolution2 * i].y, velocity0[j + resolution2 * i].x) * displayVelocityMultiplier * 0.1f;
                    particles_vectorField[ii].startColor = new Color(1, 1, 1, 1);
                    ++ii;
                }
            particleSystem_vectorField.SetParticles(particles_vectorField, voxelCount);
        }
        else renderer_vectorField.enabled = false;
    }

    UnityEngine.Vector3 IdxToPos(int idx)
    {
        if (idx >= 0 && idx < resolution2 * resolution2) return new UnityEngine.Vector3(voxelScale * Mathf.Floor(idx / resolution), (voxelScale * (idx - resolution * Mathf.Floor(idx / resolution))), 0);//new UnityEngine.Vector3(voxelScale * (idx - resolution * Mathf.Floor(idx / resolution)), voxelScale * Mathf.Floor(idx / resolution), 0);
        else return new UnityEngine.Vector3(-1, -1, -1);
    }

    int PosToIdx(UnityEngine.Vector3 pos)
    {
        UnityEngine.Vector3 scale = 2f * transform.localScale;
        if (pos.x >= 0 && pos.x < scale.x / transform.localScale.x && pos.y >= 0 && pos.y < scale.y / transform.localScale.y)
        {
            return (int)Mathf.Floor(pos.x * resolution) + (int)(Mathf.Floor(pos.y * resolution) * resolution);
        }
        else return -1;
    }

    Complex J(Complex x, float y)
    {
        return x + y * y / x;
    }

    void Joukowsky()
    {
        if (airfoilType == DisplayGrid.joukowsky)
        {
            if (xDisplacement0 != xDisplacement || yDisplacement0 != yDisplacement || angleOfAttack0 != angleOfAttack || (int)displayJoukowsky0 != (int)displayJoukowsky || airfoilType0 != airfoilType)
            {
                int ii = 0;
                float a = resolution * resolution * Time.deltaTime;
                angleOfAttack0 = angleOfAttack;
                xDisplacement0 = xDisplacement;
                yDisplacement0 = yDisplacement;
                angleOfAttack0 = angleOfAttack;
                displayJoukowsky0 = displayJoukowsky;
                airfoilType0 = airfoilType;
                float alpha = -angleOfAttack * Mathf.PI / 180f;
                Complex s = new Complex(-xDisplacement, yDisplacement);
                float k = 2f * (float)Complex.Abs(s - 1) * freestreamVelocity * Mathf.Sin(alpha);
                float circulation = k / (2f * Mathf.PI);
                float lambda = 1 + (float)s.Magnitude;
                Complex flow, zeta;
                double min = -1.65, max = 1.65;

                for (double i = min + 1.3; i <= max + 1.3 + (Math.Abs(min) + Math.Abs(max)) / resolution; i += (Math.Abs(min) + Math.Abs(max)) / resolution)
                {
                    for (double j = min; j <= max + (Math.Abs(min) + Math.Abs(max)) / resolution; j += (Math.Abs(min) + Math.Abs(max)) / resolution)
                    {
                        Complex z = 1.5f * new Complex(j, i - 0.2);
                        if (j < 0) zeta = 0.5 * (z - Complex.Sqrt(z * z - 4));
                        else if (j > 0) zeta = 0.5 * (z + Complex.Sqrt(z * z - 4));
                        flow = (freestreamVelocity * Complex.Exp(Complex.ImaginaryOne * alpha)) - (freestreamVelocity * Complex.Pow(Complex.Abs(s - 1), 2) * Complex.Exp(-Complex.ImaginaryOne * alpha) / ((zeta - s) * (zeta - s))) + (Complex.ImaginaryOne * -k) / (zeta - s);
                        joukowskyPressure[ii] = 1 - ((float)flow.Real * (float)flow.Real) / (freestreamVelocity * freestreamVelocity);
                        if (displayJoukowsky == DisplayJoukowsky.analyticVelocity) velocityMagnitude[ii] = 0.75f * (float)flow.Real;
                        else if (displayJoukowsky == DisplayJoukowsky.analyticPressure) joukowskyPressure[ii] = 1 - ((float)flow.Real * (float)flow.Real) / (freestreamVelocity * freestreamVelocity);
                        ++ii;
                    }
                }

                for (int i = 0; i <= vxlcnt; ++i)
                {
                    if (i < vxlcnt)
                    {
                        Complex circle = Complex.Exp(Complex.ImaginaryOne * Math.PI * 0.02 * i) + s;
                        Complex airfoilTransform = J(circle, lambda) * new Complex(Mathf.Cos(alpha), -Mathf.Sin(alpha));
                        if (displayJoukowsky == DisplayJoukowsky.analyticVelocity || displayJoukowsky == DisplayJoukowsky.analyticPressure) vertices[i] = 0.9f * new UnityEngine.Vector3(-(chord * (0.325f - ((xDisplacement + yDisplacement) * 0.1f))) * (float)airfoilTransform.Real, 0.75f * chord * (0.25f - ((xDisplacement + yDisplacement) * 0.1f)) * (float)airfoilTransform.Imaginary, 0);
                        else if (displayJoukowsky == DisplayJoukowsky.numericSolver) vertices[i] = 0.5f * new UnityEngine.Vector3(-(chord * (0.325f - ((xDisplacement + yDisplacement) * 0.1f))) * (float)airfoilTransform.Real, chord * (0.25f - ((xDisplacement + yDisplacement) * 0.1f)) * (float)airfoilTransform.Imaginary, 0);
                    }

                    else if (i == vxlcnt) vertices[i] = new UnityEngine.Vector3(0.5f * xDisplacement, 0);
                }

                airfoil.transform.localPosition = new UnityEngine.Vector3(1f, 1f, 0);
                boxCollider.center = new UnityEngine.Vector3(-0.01f, -0.02f, 0);
                airfoil.transform.GetChild(0).GetChild(0).position = new UnityEngine.Vector3(0, 0, 0.01f);
                airfoil.transform.GetChild(0).GetChild(1).position = new UnityEngine.Vector3(0, 0.9f, 0);
                airfoil.mesh.vertices = vertices;
            }
            else return;
        }
    }

    void NacaFourDigit()
    {
        if (airfoilType == DisplayGrid.nacaFourDigit)
        {
            if (maxCamber40 != maxCamber4 || maxCamberPosition40 != maxCamberPosition4 || maxThickness40 != maxThickness4 || angleOfAttack0 != angleOfAttack || airfoilType0 != airfoilType)
            {
                maxCamber40 = maxCamber4;
                maxCamberPosition40 = maxCamberPosition4;
                maxThickness40 = maxThickness4;
                angleOfAttack0 = angleOfAttack;
                maxCamberActual = maxCamber4 * chordActual * 0.01f;
                maxCamberPositionActual = maxCamberPosition4 * chordActual * 0.1f;
                thicknessActual = maxThickness4 * chordActual * 0.01f;
                airfoilType0 = airfoilType;
                float alpha = -angleOfAttack * Mathf.PI / 180f;

                for (int i = 0; i <= vxlcnt; ++i)
                {
                    if (maxCamber4 != 0 || maxCamberPosition4 != 0)
                    {
                        if (position >= 0 && position <= maxCamberPositionActual * chordActual)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }
                        }

                        else if (position > maxCamberPositionActual * chordActual && position <= chordActual)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * ((1 - 2 * maxCamberPositionActual) + 2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * ((1 - 2 * maxCamberPositionActual) + 2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }
                        }
                    }

                    else camberFunction = derivativeCamberFunction = 0;

                    angle = Mathf.Atan(derivativeCamberFunction);

                    if (i < hlfvxlcnt)
                    {
                        position = chordActual - (chordActual * i / hlfvxlcnt);
                        thickness = (float)(5 * thicknessActual * ((0.2969 * Mathf.Sqrt(position)) - (0.1260 * position) - (0.3516 * Mathf.Pow(position, 2)) + (0.2843 * Mathf.Pow(position, 3)) - (0.1036 * Mathf.Pow(position, 4))));
                        vertices[i] = 0.75f * new UnityEngine.Vector3((position - thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (camberFunction + thickness * Mathf.Cos(angle)) * Mathf.Sin(alpha) - 0.2f, (position - thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (camberFunction + thickness * Mathf.Cos(angle)) * Mathf.Cos(alpha) - 0.5f * alpha - 0.7f * maxCamberActual);
                    }

                    else if (i >= hlfvxlcnt && i < vxlcnt)
                    {
                        position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual)) + 0.01f;
                        thickness = (float)(5 * thicknessActual * ((0.2969 * Mathf.Sqrt(position)) - (0.1260 * position) - (0.3516 * Mathf.Pow(position, 2)) + (0.2843 * Mathf.Pow(position, 3)) - (0.1036 * Mathf.Pow(position, 4))));
                        vertices[i] = 0.75f * new UnityEngine.Vector3((position + thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (camberFunction - thickness * Mathf.Cos(angle)) * Mathf.Sin(alpha) - 0.2f, (position + thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (camberFunction - thickness * Mathf.Cos(angle)) * Mathf.Cos(alpha) - 0.5f * alpha - 0.7f * maxCamberActual);
                    }

                    else if (i == vxlcnt) vertices[i] = new UnityEngine.Vector3(0.75f * maxCamberPositionActual * Mathf.Cos(alpha) - 0.75f * maxCamberActual * Mathf.Sin(alpha) - 0.2f, 0.75f * maxCamberPositionActual * Mathf.Sin(alpha) + 0.75f * maxCamberActual * Mathf.Cos(alpha) - 0.4f * alpha - 0.5f * maxCamberActual);
                }

                airfoil.transform.position = new UnityEngine.Vector3(0.75f, 1f, 0);
                airfoil.transform.GetChild(0).GetChild(0).position = new UnityEngine.Vector3(0, 0, 0.01f);
                airfoil.transform.GetChild(0).GetChild(1).position = new UnityEngine.Vector3(0, 0.9f, 0);
                boxCollider.center = new UnityEngine.Vector3(0.2375f, -0.02f, 0);
                airfoil.mesh.vertices = vertices;
            }
            else return;
        }
    }

    void NacaFourDigitModified()
    {
        if (airfoilType == DisplayGrid.nacaFourDigitModified)
        {
            if (maxCamber4m0 != maxCamber4m || maxCamberPosition4m0 != maxCamberPosition4m || maxThickness4m0 != maxThickness4m || angleOfAttack0 != angleOfAttack || modifications4m0 != modifications4m || airfoilType0 != airfoilType)
            {
                maxCamber4m0 = maxCamber4m;
                maxCamberPosition4m0 = maxCamberPosition4m;
                maxThickness4m0 = maxThickness4m;
                angleOfAttack0 = angleOfAttack;
                maxCamberActual = maxCamber4m * chord * 0.01f;
                maxCamberPositionActual = maxCamberPosition4m * chord * 0.1f;
                thicknessActual = maxThickness4m * chord * 0.01f;
                modifications4m0 = modifications4m;
                airfoilType0 = airfoilType;
                float alpha = -angleOfAttack * Mathf.PI / 180f;

                for (int i = 0; i <= vxlcnt; ++i)
                {
                    if (maxCamber4m != 0 || maxCamberPosition4m != 0)
                    {
                        if (position >= 0 && position <= maxCamberPositionActual * chordActual)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }
                        }

                        else if (position > maxCamberPositionActual * chord && position <= chordActual)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * ((1 - 2 * maxCamberPositionActual) + 2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * ((1 - 2 * maxCamberPositionActual) + 2 * maxCamberPositionActual * (position / chordActual) - Mathf.Pow(position / chordActual, 2));
                                derivativeCamberFunction = 2 * maxCamberActual / Mathf.Pow(1 - maxCamberPositionActual, 2) * (maxCamberPositionActual - (position / chordActual));
                            }
                        }
                    }

                    else camberFunction = derivativeCamberFunction = 0;

                    angle = Mathf.Atan(derivativeCamberFunction);

                    if (position >= 0 && position <= thicknessActual)
                    {
                        if (modifications4m == Display4m.Three)
                        {
                            a0 = 0;
                            a1 = 0.920286f;
                            a2 = -2.801900f;
                            a3 = 2.817990f;
                        }

                        else if (modifications4m == Display4m.Five)
                        {
                            a0 = 0;
                            a1 = 0.477f;
                            a2 = -0.708f;
                            a3 = 0.308f;
                        }

                        else if (modifications4m == Display4m.ThirtyThree)
                        {
                            a0 = 0.14845f;
                            a1 = 0.412103f;
                            a2 = -1.67261f;
                            a3 = 1.68869f;
                        }

                        else if (modifications4m == Display4m.ThirtyFour)
                        {
                            a0 = 0.14845f;
                            a1 = 0.193233f;
                            a2 = -0.558166f;
                            a3 = 0.283208f;
                        }

                        else if (modifications4m == Display4m.ThirtyFive)
                        {
                            a0 = 0.14845f;
                            a1 = 0.083362f;
                            a2 = -0.18315f;
                            a3 = -0.00691f;
                        }

                        else if (modifications4m == Display4m.SixtyTwo)
                        {
                            a0 = 0.2969f;
                            a1 = 0.213337f;
                            a2 = -2.931954f;
                            a3 = 5.22917f;
                        }

                        else if (modifications4m == Display4m.SixtyThree)
                        {
                            a0 = 0.2969f;
                            a1 = -0.096082f;
                            a2 = -0.543310f;
                            a3 = 0.559395f;
                        }

                        else if (modifications4m == Display4m.SixtyFour)
                        {
                            a0 = 0.2969f;
                            a1 = -0.246867f;
                            a2 = 0.175384f;
                            a3 = -0.266917f;
                        }

                        else if (modifications4m == Display4m.SixtyFive)
                        {
                            a0 = 0.2969f;
                            a1 = -0.310275f;
                            a2 = 0.3417f;
                            a3 = -0.321820f;
                        }

                        else if (modifications4m == Display4m.NinetyThree)
                        {
                            a0 = 0.514246f;
                            a1 = -0.840115f;
                            a2 = 1.1101f;
                            a3 = -1.09401f;
                        }

                        if (i < hlfvxlcnt)
                        {
                            position = chordActual - (chordActual * i / hlfvxlcnt);
                            thickness = a0 * Mathf.Sqrt(position) + a1 * position + a2 * Mathf.Pow(position, 2) + a3 * Mathf.Pow(position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position - thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position - thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha - 0.5f * thicknessActual + 0.05f);
                        }

                        else if (i >= hlfvxlcnt && i < vxlcnt)
                        {
                            position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                            thickness = a0 * Mathf.Sqrt(position) + a1 * position + a2 * Mathf.Pow(position, 2) + a3 * Mathf.Pow(position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position + thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position + thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha - 0.5f * thicknessActual + 0.05f);
                        }
                    }

                    else if (position > thicknessActual && position <= 1)
                    {
                        if (modifications4m == Display4m.Three)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        else if (modifications4m == Display4m.Five)
                        {
                            d0 = 0.002f;
                            d1 = 0.465f;
                            d2 = -0.684f;
                            d3 = 0.292f;
                        }

                        else if (modifications4m == Display4m.ThirtyThree)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        else if (modifications4m == Display4m.ThirtyFour)
                        {
                            d0 = 0.002f;
                            d1 = 0.315f;
                            d2 = -0.233333f;
                            d3 = -0.032407f;
                        }

                        else if (modifications4m == Display4m.ThirtyFive)
                        {
                            d0 = 0.002f;
                            d1 = 0.465f;
                            d2 = -0.684f;
                            d3 = 0.292f;
                        }

                        else if (modifications4m == Display4m.SixtyTwo)
                        {
                            d0 = 0.002f;
                            d1 = 0.2f;
                            d2 = -0.040625f;
                            d3 = -0.070312f;
                        }

                        else if (modifications4m == Display4m.SixtyThree)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        else if (modifications4m == Display4m.SixtyFour)
                        {
                            d0 = 0.002f;
                            d1 = 0.315f;
                            d2 = -0.233333f;
                            d3 = -0.032407f;
                        }

                        else if (modifications4m == Display4m.SixtyFive)
                        {
                            d0 = 0.002f;
                            d1 = 0.465f;
                            d2 = -0.684f;
                            d3 = 0.292f;
                        }

                        else if (modifications4m == Display4m.NinetyThree)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        if (i < hlfvxlcnt)
                        {
                            position = chordActual - (chordActual * i / hlfvxlcnt);
                            thickness = d0 + d1 * (1 - position) + d2 * Mathf.Pow(1 - position, 2) + d3 * Mathf.Pow(1 - position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position - thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position - thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha - 0.5f * thicknessActual + 0.05f);
                        }

                        else if (i >= hlfvxlcnt && i < vxlcnt)
                        {
                            position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                            thickness = d0 + d1 * (1 - position) + d2 * Mathf.Pow(1 - position, 2) + d3 * Mathf.Pow(1 - position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position + thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position + thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha - 0.5f * thicknessActual + 0.05f);
                        }
                    }

                    if (i == vxlcnt)
                    {
                        float x = 1.5f * 0.5f * maxCamberPositionActual;
                        float y = 1.5f * 1.75f * maxCamberActual * thicknessActual + 0.2f * thicknessActual - 0.025f * maxCamberPositionActual;
                        vertices[i] = new UnityEngine.Vector3(x * Mathf.Cos(alpha) - y * Mathf.Sin(alpha) - 0.2f, x * Mathf.Sin(alpha) + y * Mathf.Cos(alpha) - 0.4f * alpha - 0.5f * thicknessActual + 0.05f);
                    }
                }

                airfoil.transform.position = new UnityEngine.Vector3(0.75f, 1f, 0);
                airfoil.transform.GetChild(0).GetChild(0).position = new UnityEngine.Vector3(0, 0, 0.01f);
                airfoil.transform.GetChild(0).GetChild(1).position = new UnityEngine.Vector3(0, 0.9f, 0);
                boxCollider.center = new UnityEngine.Vector3(0.2375f, -0.02f, 0);
                airfoil.mesh.vertices = vertices;
            }

            else return;
        }
    }

    void NacaFiveDigit()
    {
        if (airfoilType == DisplayGrid.nacaFiveDigit)
        {
            if (designLiftCoefficient50 != designLiftCoefficient5 || maxCamberPosition50 != maxCamberPosition5 || camberType50 != camberType5 || maxThickness50 != maxThickness5 || angleOfAttack0 != angleOfAttack || airfoilType0 != airfoilType)
            {
                designLiftCoefficient50 = designLiftCoefficient5;
                maxCamberPosition50 = maxCamberPosition5;
                camberType50 = camberType5;
                maxThickness50 = maxThickness5;
                angleOfAttack0 = angleOfAttack;
                designLiftCoefficientActual = designLiftCoefficient5 * 0.15f;
                maxCamberPositionActual = maxCamberPosition5 * 0.05f * chordActual;
                thicknessActual = maxThickness5 * 0.01f * chordActual;
                airfoilType0 = airfoilType;
                float alpha = -angleOfAttack * Mathf.PI / 180f;

                for (int i = 0; i <= vxlcnt; ++i)
                {
                    if (camberType5 == 0)
                    {
                        if (maxCamberPosition5 == 1)
                        {
                            r = 0.058f;
                            k1 = 361.4f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5 == 2)
                        {
                            r = 0.1260f;
                            k1 = 51.640f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5 == 3)
                        {
                            r = 0.2025f;
                            k1 = 15.957f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5 == 4)
                        {
                            r = 0.2900f;
                            k1 = 6.643f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5 == 5)
                        {
                            r = 0.3910f;
                            k1 = 3.230f * designLiftCoefficientActual / 0.3f;
                        }

                        if (position >= 0 && position < r)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 / 6f * (Mathf.Pow(position, 3) - 3 * r * Mathf.Pow(position, 2) + r * r * (3 - r) * position);
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position, 2) - 6 * r * position + r * r * (3 - r));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 / 6f * (Mathf.Pow(position, 3) - 3 * r * Mathf.Pow(position, 2) + r * r * (3 - r) * position);
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position, 2) - 6 * r * position + r * r * (3 - r));
                            }
                        }

                        else if (position >= r && position < 1)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 * Mathf.Pow(r, 3) / 6f * (1 - position);
                                derivativeCamberFunction = -k1 * Mathf.Pow(r, 3) / 6f;
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 * Mathf.Pow(r, 3) / 6f * (1 - position);
                                derivativeCamberFunction = -k1 * Mathf.Pow(r, 3) / 6f;
                            }
                        }
                    }

                    else if (camberType5 == 1)
                    {
                        if (maxCamberPosition5 == 1) maxCamberPosition5 = 2;

                        if (maxCamberPosition5 == 2)
                        {
                            r = 0.1300f;
                            k1 = 51.990f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.000764f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5 == 3)
                        {
                            r = 0.2170f;
                            k1 = 15.793f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.00677f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5 == 4)
                        {
                            r = 0.3180f;
                            k1 = 6.520f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.0303f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5 == 5)
                        {
                            r = 0.4410f;
                            k1 = 3.191f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.1355f * designLiftCoefficientActual / 0.3f;
                        }

                        if (position >= 0 && position < r)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 / 6f * (Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 / 6f * (Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }
                        }

                        else if (position >= r && position <= 1)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 / 6f * (k1k2 * Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3 * k1k2 * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 / 6f * (k1k2 * Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3 * k1k2 * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }
                        }
                    }

                    angle = Mathf.Atan(derivativeCamberFunction);

                    if (i < hlfvxlcnt)
                    {
                        position = chordActual - (chordActual * i / hlfvxlcnt);
                        thickness = (float)(5 * thicknessActual * ((0.2969 * Mathf.Sqrt(position)) - (0.1260 * position) - (0.3516 * Mathf.Pow(position, 2)) + (0.2843 * Mathf.Pow(position, 3)) - (0.1036 * Mathf.Pow(position, 4))));
                        vertices[i] = 0.75f * new UnityEngine.Vector3((position - thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (camberFunction + thickness * Mathf.Cos(angle)) * Mathf.Sin(alpha) - 0.2f, (position - thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (camberFunction + thickness * Mathf.Cos(angle)) * Mathf.Cos(alpha) - 0.5f * alpha);
                    }

                    else if (i >= hlfvxlcnt && i < vxlcnt)
                    {
                        position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                        thickness = (float)(5 * thicknessActual * ((0.2969 * Mathf.Sqrt(position)) - (0.1260 * position) - (0.3516 * Mathf.Pow(position, 2)) + (0.2843 * Mathf.Pow(position, 3)) - (0.1036 * Mathf.Pow(position, 4))));
                        vertices[i] = 0.75f * new UnityEngine.Vector3((position + thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (camberFunction - thickness * Mathf.Cos(angle)) * Mathf.Sin(alpha) - 0.2f, (position + thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (camberFunction - thickness * Mathf.Cos(angle)) * Mathf.Cos(alpha) - 0.5f * alpha);
                    }

                    else if (i == vxlcnt) vertices[i] = new UnityEngine.Vector3(0.75f * maxCamberPositionActual * Mathf.Cos(alpha) - 0.2f, 0.75f * maxCamberPositionActual * Mathf.Sin(alpha) - 0.4f * alpha);
                }

                airfoil.transform.position = new UnityEngine.Vector3(0.75f, 1f, 0);
                airfoil.transform.GetChild(0).GetChild(0).position = new UnityEngine.Vector3(0, 0, 0.01f);
                airfoil.transform.GetChild(0).GetChild(1).position = new UnityEngine.Vector3(0, 0.9f, 0);
                boxCollider.center = new UnityEngine.Vector3(0.2375f, -0.02f, 0);
                airfoil.mesh.vertices = vertices;
            }

            else return;
        }
    }

    void NacaFiveDigitModified()
    {
        if (airfoilType == DisplayGrid.nacaFiveDigitModified)
        {
            if (designLiftCoefficient5m0 != designLiftCoefficient5m || maxCamberPosition5m0 != maxCamberPosition5m || camberType5m0 != camberType5m || maxThickness5m0 != maxThickness5m || angleOfAttack0 != angleOfAttack || modifications5m0 != modifications5m || airfoilType0 != airfoilType)
            {
                designLiftCoefficient5m0 = designLiftCoefficient5m;
                maxCamberPosition5m0 = maxCamberPosition5m;
                camberType5m0 = camberType5m;
                maxThickness5m0 = maxThickness5m;
                angleOfAttack0 = angleOfAttack;
                designLiftCoefficientActual = designLiftCoefficient5m * 0.15f;
                maxCamberPositionActual = maxCamberPosition5m * 0.05f * chordActual;
                thicknessActual = maxThickness5m * 0.01f * chordActual;
                modifications5m0 = modifications5m;
                airfoilType0 = airfoilType;
                float alpha = -angleOfAttack * Mathf.PI / 180f;

                for (int i = 0; i <= vxlcnt; ++i)
                {
                    if (camberType5m == 0)
                    {
                        if (maxCamberPosition5m == 1)
                        {
                            r = 0.058f;
                            k1 = 361.4f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5m == 2)
                        {
                            r = 0.1260f;
                            k1 = 51.640f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5m == 3)
                        {
                            r = 0.2025f;
                            k1 = 15.957f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5m == 4)
                        {
                            r = 0.2900f;
                            k1 = 6.643f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5m == 5)
                        {
                            r = 0.3910f;
                            k1 = 3.230f * designLiftCoefficientActual / 0.3f;
                        }

                        if (position >= 0 && position < r)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 / 6f * (Mathf.Pow(position, 3) - 3 * r * Mathf.Pow(position, 2) + r * r * (3 - r) * position);
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position, 2) - 6 * r * position + r * r * (3 - r));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 / 6f * (Mathf.Pow(position, 3) - 3 * r * Mathf.Pow(position, 2) + r * r * (3 - r) * position);
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position, 2) - 6 * r * position + r * r * (3 - r));
                            }
                        }

                        else if (position >= r && position < 1)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 * Mathf.Pow(r, 3) / 6f * (1 - position);
                                derivativeCamberFunction = -k1 * Mathf.Pow(r, 3) / 6f;
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 * Mathf.Pow(r, 3) / 6f * (1 - position);
                                derivativeCamberFunction = -k1 * Mathf.Pow(r, 3) / 6f;
                            }
                        }
                    }

                    else if (camberType5m == 1)
                    {
                        if (maxCamberPosition5m == 1) maxCamberPosition5m = 2;

                        if (maxCamberPosition5m == 2)
                        {
                            r = 0.1300f;
                            k1 = 51.990f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.000764f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5m == 3)
                        {
                            r = 0.2170f;
                            k1 = 15.793f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.00677f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5m == 4)
                        {
                            r = 0.3180f;
                            k1 = 6.520f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.0303f * designLiftCoefficientActual / 0.3f;
                        }

                        else if (maxCamberPosition5m == 5)
                        {
                            r = 0.4410f;
                            k1 = 3.191f * designLiftCoefficientActual / 0.3f;
                            k1k2 = 0.1355f * designLiftCoefficientActual / 0.3f;
                        }

                        if (position >= 0 && position < r)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 / 6f * (Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 / 6f * (Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3f * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }
                        }

                        else if (position >= r && position <= 1)
                        {
                            if (i < hlfvxlcnt)
                            {
                                position = chordActual - (chordActual * i / hlfvxlcnt);
                                camberFunction = k1 / 6f * (k1k2 * Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3 * k1k2 * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }

                            else if (i >= hlfvxlcnt && i < vxlcnt)
                            {
                                position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                                camberFunction = k1 / 6f * (k1k2 * Mathf.Pow(position - r, 3) - k1k2 * Mathf.Pow(1 - r, 3) * position - Mathf.Pow(r, 3) * position + Mathf.Pow(r, 3));
                                derivativeCamberFunction = k1 / 6f * (3 * k1k2 * Mathf.Pow(position - r, 2) - k1k2 * Mathf.Pow(1 - r, 3) - Mathf.Pow(r, 3));
                            }
                        }
                    }

                    angle = Mathf.Atan(derivativeCamberFunction);

                    if (position >= 0 && position <= thicknessActual)
                    {
                        if (modifications5m == Display5m.Three)
                        {
                            a0 = 0;
                            a1 = 0.920286f;
                            a2 = -2.801900f;
                            a3 = 2.817990f;
                        }

                        else if (modifications5m == Display5m.Five)
                        {
                            a0 = 0;
                            a1 = 0.477f;
                            a2 = -0.708f;
                            a3 = 0.308f;
                        }

                        else if (modifications5m == Display5m.ThirtyThree)
                        {
                            a0 = 0.14845f;
                            a1 = 0.412103f;
                            a2 = -1.67261f;
                            a3 = 1.68869f;
                        }

                        else if (modifications5m == Display5m.ThirtyFour)
                        {
                            a0 = 0.14845f;
                            a1 = 0.193233f;
                            a2 = -0.558166f;
                            a3 = 0.283208f;
                        }

                        else if (modifications5m == Display5m.ThirtyFive)
                        {
                            a0 = 0.14845f;
                            a1 = 0.083362f;
                            a2 = -0.18315f;
                            a3 = -0.00691f;
                        }

                        else if (modifications5m == Display5m.SixtyTwo)
                        {
                            a0 = 0.2969f;
                            a1 = 0.213337f;
                            a2 = -2.931954f;
                            a3 = 5.22917f;
                        }

                        else if (modifications5m == Display5m.SixtyThree)
                        {
                            a0 = 0.2969f;
                            a1 = -0.096082f;
                            a2 = -0.543310f;
                            a3 = 0.559395f;
                        }

                        else if (modifications5m == Display5m.SixtyFour)
                        {
                            a0 = 0.2969f;
                            a1 = -0.246867f;
                            a2 = 0.175384f;
                            a3 = -0.266917f;
                        }

                        else if (modifications5m == Display5m.SixtyFive)
                        {
                            a0 = 0.2969f;
                            a1 = -0.310275f;
                            a2 = 0.3417f;
                            a3 = -0.321820f;
                        }

                        else if (modifications5m == Display5m.NinetyThree)
                        {
                            a0 = 0.514246f;
                            a1 = -0.840115f;
                            a2 = 1.1101f;
                            a3 = -1.09401f;
                        }

                        if (i < hlfvxlcnt)
                        {
                            position = chordActual - (chordActual * i / hlfvxlcnt);
                            thickness = a0 * Mathf.Sqrt(position) + a1 * position + a2 * Mathf.Pow(position, 2) + a3 * Mathf.Pow(position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position - thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position - thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha);
                        }

                        else if (i >= hlfvxlcnt && i < vxlcnt)
                        {
                            position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                            thickness = a0 * Mathf.Sqrt(position) + a1 * position + a2 * Mathf.Pow(position, 2) + a3 * Mathf.Pow(position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position + thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position + thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha);
                        }
                    }

                    else if (position > thicknessActual && position <= 1)
                    {
                        if (modifications5m == Display5m.Three)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        else if (modifications5m == Display5m.Five)
                        {
                            d0 = 0.002f;
                            d1 = 0.465f;
                            d2 = -0.684f;
                            d3 = 0.292f;
                        }

                        else if (modifications5m == Display5m.ThirtyThree)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        else if (modifications5m == Display5m.ThirtyFour)
                        {
                            d0 = 0.002f;
                            d1 = 0.315f;
                            d2 = -0.233333f;
                            d3 = -0.032407f;
                        }

                        else if (modifications5m == Display5m.ThirtyFive)
                        {
                            d0 = 0.002f;
                            d1 = 0.465f;
                            d2 = -0.684f;
                            d3 = 0.292f;
                        }

                        else if (modifications5m == Display5m.SixtyTwo)
                        {
                            d0 = 0.002f;
                            d1 = 0.2f;
                            d2 = -0.040625f;
                            d3 = -0.070312f;
                        }

                        else if (modifications5m == Display5m.SixtyThree)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        else if (modifications5m == Display5m.SixtyFour)
                        {
                            d0 = 0.002f;
                            d1 = 0.315f;
                            d2 = -0.233333f;
                            d3 = -0.032407f;
                        }

                        else if (modifications5m == Display5m.SixtyFive)
                        {
                            d0 = 0.002f;
                            d1 = 0.465f;
                            d2 = -0.684f;
                            d3 = 0.292f;
                        }

                        else if (modifications5m == Display5m.NinetyThree)
                        {
                            d0 = 0.002f;
                            d1 = 0.234f;
                            d2 = -0.068571f;
                            d3 = -0.093878f;
                        }

                        if (i < hlfvxlcnt)
                        {
                            position = chordActual - (chordActual * i / hlfvxlcnt);
                            thickness = d0 + d1 * (1 - position) + d2 * Mathf.Pow(1 - position, 2) + d3 * Mathf.Pow(1 - position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position - thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position - thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction + thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha);
                        }

                        else if (i >= hlfvxlcnt && i < vxlcnt)
                        {
                            position = 2f * ((chordActual * i / vxlcnt) - (0.5f * chordActual));
                            thickness = d0 + d1 * (1 - position) + d2 * Mathf.Pow(1 - position, 2) + d3 * Mathf.Pow(1 - position, 3);
                            vertices[i] = 0.75f * new UnityEngine.Vector3((position + thickness * Mathf.Sin(angle)) * Mathf.Cos(alpha) - (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Sin(alpha) - 0.2f, (position + thickness * Mathf.Sin(angle)) * Mathf.Sin(alpha) + (thicknessActual * 5f * (camberFunction - thickness * Mathf.Cos(angle))) * Mathf.Cos(alpha) - 0.5f * alpha);
                        }
                    }

                    if (i == vxlcnt) vertices[i] = new UnityEngine.Vector3(0.75f * maxCamberPositionActual * Mathf.Cos(alpha) - 0.2f, 0.75f * maxCamberPositionActual * Mathf.Sin(alpha) - 0.4f * alpha);
                }

                airfoil.transform.position = new UnityEngine.Vector3(0.75f, 1f, 0);
                airfoil.transform.GetChild(0).GetChild(0).position = new UnityEngine.Vector3(0, 0, 0.01f);
                airfoil.transform.GetChild(0).GetChild(1).position = new UnityEngine.Vector3(0, 0.9f, 0);
                boxCollider.center = new UnityEngine.Vector3(0.2375f, -0.02f, 0);
                airfoil.mesh.vertices = vertices;
            }
            else return;
        }
    }

    int vxlcnt, hlfvxlcnt;

    void Start()
    {
        airfoil = GetComponent<MeshFilter>();
        chordActual = 1;
        vertices = airfoil.sharedMesh.vertices;
        vxlcnt = vertices.Length - 1;
        hlfvxlcnt = vxlcnt / 2;
        circle = airfoil.sharedMesh.vertices;
        boxCollider = GetComponent<BoxCollider>();
        particleSystem_scalarField = transform.GetChild(0).GetChild(0).GetComponent<ParticleSystem>();
        renderer_scalarField = particleSystem_scalarField.GetComponent<Renderer>();
        particleSystem_vectorField = transform.GetChild(0).GetChild(1).GetComponent<ParticleSystem>();
        renderer_vectorField = particleSystem_vectorField.GetComponent<Renderer>();
    }

    void Update()
    {
        Initialize();
        Joukowsky();
        NacaFourDigit();
        NacaFourDigitModified();
        NacaFiveDigit();
        NacaFiveDigitModified();
        Emit();
        Diffuse();
        Project();
        AdvectVelocity();
        Project();
        AdvectDensity();
        VorticityConfinement();
        Display();
    }
}
