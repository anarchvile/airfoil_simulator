using System;
using UnityEngine;
using UnityEngine.UI;

public class ui : MonoBehaviour
{
    public solver solver;
    [Header("ranges")]
    public Vector2 resolution = new Vector2(10, 400);
    public Vector2 inletRadius = new Vector2(0, 1);
    public Vector2 inletSamples = new Vector2(0, 100);
    public Vector2 radius = new Vector2(0, 1);
    public Vector2 samples = new Vector2(0, 100);
    public Vector2 friction = new Vector2(0, 1);
    public Vector2 acceleration = new Vector2(-1, 1);
    public Vector2 swirl = new Vector2(0, 1);
    public Vector2 densityDiffusion = new Vector2(0, 1);
    public Vector2 velocityDiffusion = new Vector2(0, 1);
    public Vector2 diffusionQuality = new Vector2(0, 10);
    public Vector2 solverQuality = new Vector2(1, 50);
    public Vector2 angleOfAttack = new Vector2(0, 20);
    public Vector2 xDisplacement = new Vector2(0, 0.25f);
    public Vector2 yDisplacement = new Vector2(0, 0.25f);
    public Vector2 maxCamberPosition4 = new Vector2(0, 10);

    public void Start()
    {
        //framerate = transform.GetChild(0).GetChild(transform.GetChild(0).childCount - 1).GetChild(0).GetComponent<Text>();

        Slider slider;

        transform.GetChild(0).GetChild(0).GetChild(0).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().resolution.ToString();
        slider = transform.GetChild(0).GetChild(0).GetChild(0).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().resolution;
        slider.minValue = resolution.x;
        slider.maxValue = resolution.y;

        transform.GetChild(0).GetChild(0).GetChild(1).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().acceleration.ToString();
        slider = transform.GetChild(0).GetChild(0).GetChild(1).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().acceleration;
        slider.minValue = 0;
        slider.maxValue = 10;
        
        transform.GetChild(0).GetChild(0).GetChild(3).GetChild(0).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().injectInletDensity.ToString();

        transform.GetChild(0).GetChild(0).GetChild(3).GetChild(1).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().inletRadius.ToString();
        slider = transform.GetChild(0).GetChild(0).GetChild(3).GetChild(1).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().inletRadius;
        slider.minValue = inletRadius.x;
        slider.maxValue = inletRadius.y;

        transform.GetChild(0).GetChild(0).GetChild(3).GetChild(2).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().inletSamples.ToString();
        slider = transform.GetChild(0).GetChild(0).GetChild(3).GetChild(2).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().inletSamples;
        slider.minValue = inletSamples.x;
        slider.maxValue = inletSamples.y;

        transform.GetChild(0).GetChild(0).GetChild(4).GetChild(0).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().injectDensity.ToString();

        transform.GetChild(0).GetChild(0).GetChild(4).GetChild(1).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().injectVelocity.ToString();

        transform.GetChild(0).GetChild(0).GetChild(4).GetChild(2).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxInjectVelocity.ToString();

        transform.GetChild(0).GetChild(0).GetChild(4).GetChild(3).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().radius.ToString();
        slider = transform.GetChild(0).GetChild(0).GetChild(4).GetChild(3).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().radius;
        slider.minValue = radius.x;
        slider.maxValue = radius.y;

        transform.GetChild(0).GetChild(0).GetChild(4).GetChild(4).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().samples.ToString();
        slider = transform.GetChild(0).GetChild(0).GetChild(4).GetChild(4).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().samples;
        slider.minValue = samples.x;
        slider.maxValue = samples.y;

        transform.GetChild(0).GetChild(1).GetChild(0).GetComponent<Dropdown>().value = (int)GameObject.Find("wing_poly_geo").GetComponent<solver>().airfoilType;

        transform.GetChild(0).GetChild(1).GetChild(1).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().angleOfAttack.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(1).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().angleOfAttack;
        slider.minValue = angleOfAttack.x;
        slider.maxValue = angleOfAttack.y;

        transform.GetChild(0).GetChild(1).GetChild(2).GetComponent<Dropdown>().value = (int)GameObject.Find("wing_poly_geo").GetComponent<solver>().displayJoukowsky;

        transform.GetChild(0).GetChild(1).GetChild(3).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().xDisplacement.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(3).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().xDisplacement;
        slider.minValue = xDisplacement.x;
        slider.maxValue = xDisplacement.y;

        transform.GetChild(0).GetChild(1).GetChild(4).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().yDisplacement.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(4).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().yDisplacement;
        slider.minValue = yDisplacement.x;
        slider.maxValue = yDisplacement.y;

        transform.GetChild(0).GetChild(1).GetChild(5).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(5).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4;
        slider.minValue = 0;
        slider.maxValue = 25;

        transform.GetChild(0).GetChild(1).GetChild(6).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(6).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4;
        slider.minValue = 0.5f;
        slider.maxValue = 9.5f;

        transform.GetChild(0).GetChild(1).GetChild(7).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(7).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4;
        slider.minValue = 0;
        slider.maxValue = 30;

        transform.GetChild(0).GetChild(1).GetChild(8).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4m.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(8).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4m;
        slider.minValue = 0;
        slider.maxValue = 15;

        transform.GetChild(0).GetChild(1).GetChild(9).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4m.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(9).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4m;
        slider.minValue = 0.5f;
        slider.maxValue = 9.5f;

        transform.GetChild(0).GetChild(1).GetChild(10).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4m.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(10).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4m;
        slider.minValue = 0;
        slider.maxValue = 30;

        transform.GetChild(0).GetChild(1).GetChild(11).GetComponent<Dropdown>().value = (int)GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m;

        transform.GetChild(0).GetChild(1).GetChild(12).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(12).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5;
        slider.minValue = 0.05f;
        slider.maxValue = 6.67f;

        transform.GetChild(0).GetChild(1).GetChild(13).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(13).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5;
        slider.minValue = 1;
        slider.maxValue = 5;

        transform.GetChild(0).GetChild(1).GetChild(14).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(14).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5;
        slider.minValue = 0;
        slider.maxValue = 1;

        transform.GetChild(0).GetChild(1).GetChild(15).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(15).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5;
        slider.minValue = 0;
        slider.maxValue = 30;

        transform.GetChild(0).GetChild(1).GetChild(16).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5m.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(16).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5m;
        slider.minValue = 0.05f;
        slider.maxValue = 6.67f;

        transform.GetChild(0).GetChild(1).GetChild(17).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5m.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(17).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5m;
        slider.minValue = 1;
        slider.maxValue = 5;

        transform.GetChild(0).GetChild(1).GetChild(18).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5m.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(18).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5m;
        slider.minValue = 0;
        slider.maxValue = 1;

        transform.GetChild(0).GetChild(1).GetChild(19).GetChild(1).GetComponent<InputField>().text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5m.ToString();
        slider = transform.GetChild(0).GetChild(1).GetChild(19).GetChild(2).GetComponent<Slider>();
        slider.value = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5m;
        slider.minValue = 0;
        slider.maxValue = 30;

        GameObject.Find("inletEmission").transform.localScale = Vector3.one;
        GameObject.Find("mouseEmission").transform.localScale = Vector3.zero;
        GameObject.Find("displayJoukowsky").transform.localScale = Vector3.one;
        GameObject.Find("xDisplacement").transform.localScale = Vector3.one;
        GameObject.Find("yDisplacement").transform.localScale = Vector3.one;
        GameObject.Find("maxCamber4").transform.localScale = Vector3.zero;
        GameObject.Find("maxCamberPosition4").transform.localScale = Vector3.zero;
        GameObject.Find("maxThickness4").transform.localScale = Vector3.zero;
        GameObject.Find("maxCamber4m").transform.localScale = Vector3.zero;
        GameObject.Find("maxCamberPosition4m").transform.localScale = Vector3.zero;
        GameObject.Find("maxThickness4m").transform.localScale = Vector3.zero;
        GameObject.Find("modifications4m").transform.localScale = Vector3.zero;
        GameObject.Find("designLiftCoefficient5").transform.localScale = Vector3.zero;
        GameObject.Find("maxCamberPosition5").transform.localScale = Vector3.zero;
        GameObject.Find("camberType5").transform.localScale = Vector3.zero;
        GameObject.Find("maxThickness5").transform.localScale = Vector3.zero;
        GameObject.Find("designLiftCoefficient5m").transform.localScale = Vector3.zero;
        GameObject.Find("maxCamberPosition5m").transform.localScale = Vector3.zero;
        GameObject.Find("camberType5m").transform.localScale = Vector3.zero;
        GameObject.Find("maxThickness5m").transform.localScale = Vector3.zero;
        GameObject.Find("modifications5m").transform.localScale = Vector3.zero;
        GameObject.Find("joukowsky").transform.localScale = Vector3.one;
        GameObject.Find("nacaFourDigitAirfoil").transform.localScale = Vector3.zero;
        GameObject.Find("nacaFourDigitModifiedAirfoil").transform.localScale = Vector3.zero;
        GameObject.Find("nacaFiveDigitAirfoil").transform.localScale = Vector3.zero;
        GameObject.Find("nacaFiveDigitModifiedAirfoil").transform.localScale = Vector3.zero;
    }

    void Update()
    {
        if (Input.GetKeyDown("p")) GameObject.Find("wing_poly_geo").GetComponent<solver>().enabled = !GameObject.Find("wing_poly_geo").GetComponent<solver>().enabled;
    }


    void UpdateSlider(InputField field)
    {
        field.transform.parent.GetChild(2).GetComponent<Slider>().value = Convert.ToSingle(field.text);
    }

    void UpdateField(Slider slider)
    {
        slider.transform.parent.GetChild(1).GetComponent<InputField>().text = slider.value.ToString();
    }

    public void setResolutionFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().resolution = (int)slider.value;
        UpdateField(slider);
    }

    public void setResolutionFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().resolution = (int)Mathf.Clamp(Convert.ToInt32(field.text), resolution.x, resolution.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().resolution.ToString();
        UpdateSlider(field);
    }

    public void setInletDensity(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().injectInletDensity = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().injectInletDensity.ToString();
    }

    public void setInletRadiusFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().inletRadius = slider.value;
        UpdateField(slider);
    }

    public void setInletRadiusFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().inletRadius = Mathf.Clamp(Convert.ToSingle(field.text), inletRadius.x, inletRadius.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().inletRadius.ToString();
        UpdateSlider(field);
    }

    public void setInletSamplesFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().inletSamples = (int)slider.value;
        UpdateField(slider);
    }

    public void setInletSamplesFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().inletSamples = (int)Mathf.Clamp(Convert.ToInt32(field.text), inletSamples.x, inletSamples.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().inletSamples.ToString();
        UpdateSlider(field);
    }

    public void setDensity(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().injectDensity = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().injectDensity.ToString();
    }

    public void setVelocity(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().injectVelocity = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().injectVelocity.ToString();
    }

    public void setMaxVelocity(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxInjectVelocity = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxInjectVelocity.ToString();
    }

    public void setRadiusFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().radius = slider.value;
        UpdateField(slider);
    }

    public void setRadiusFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().radius = Mathf.Clamp(Convert.ToSingle(field.text), radius.x, radius.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().radius.ToString();
        UpdateSlider(field);
    }

    public void setSamplesFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().samples = (int)slider.value;
        UpdateField(slider);
    }

    public void setSamplesFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().samples = (int)Mathf.Clamp(Convert.ToInt32(field.text), samples.x, samples.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().samples.ToString();
        UpdateSlider(field);
    }

    public void setFrictionFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().friction = slider.value;
        UpdateField(slider);
    }

    public void setFrictionFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().friction = Mathf.Clamp(Convert.ToSingle(field.text), friction.x, friction.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().friction.ToString();
        UpdateSlider(field);
    }

    public void setAccelerationFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().acceleration = slider.value;
        UpdateField(slider);
    }

    public void setAccelerationFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().acceleration = Mathf.Clamp(Convert.ToSingle(field.text), 0, 10);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().acceleration.ToString();
        UpdateSlider(field);
    }

    public void setSwirlFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().swirl = slider.value;
        UpdateField(slider);
    }

    public void setSwirlFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().swirl = Mathf.Clamp(Convert.ToSingle(field.text), swirl.x, swirl.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().swirl.ToString();
        UpdateSlider(field);
    }

    public void diffuseDensity(Toggle toggle)
    {
        if (GameObject.Find("wing_poly_geo").GetComponent<solver>().densityDissipation == 0) GameObject.Find("wing_poly_geo").GetComponent<solver>().diffuseDensity = toggle.isOn;
        else GameObject.Find("wing_poly_geo").GetComponent<solver>().diffuseDensity = toggle.isOn = true;
    }

    public void diffuseVelocity(Toggle toggle)
    {
        if (GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityDissipation == 0) GameObject.Find("wing_poly_geo").GetComponent<solver>().diffuseVelocity = toggle.isOn;
        else GameObject.Find("wing_poly_geo").GetComponent<solver>().diffuseVelocity = toggle.isOn = true;
    }

    public void setDensityDiffusionFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().densityDiffusion = slider.value;
        UpdateField(slider);
    }

    public void setDensityDiffusionFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().densityDiffusion = Mathf.Clamp(Convert.ToSingle(field.text), densityDiffusion.x, densityDiffusion.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().densityDiffusion.ToString();
        UpdateSlider(field);
    }

    public void setVelocityDiffusionFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityDiffusion = slider.value;
        UpdateField(slider);
    }

    public void setVelocityDiffusionFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityDiffusion = Mathf.Clamp(Convert.ToSingle(field.text), velocityDiffusion.x, velocityDiffusion.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityDiffusion.ToString();
        UpdateSlider(field);
    }

    public void setDiffusionQualityFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().diffusionQuality = (int)slider.value;
        UpdateField(slider);
    }

    public void setDiffusionQualityFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().diffusionQuality = (int)Mathf.Clamp(Convert.ToInt32(field.text), diffusionQuality.x, diffusionQuality.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().diffusionQuality.ToString();
        UpdateSlider(field);
    }

    public void advectDensity(Toggle toggle)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().advectDensity = toggle.isOn;
    }

    public void setVelocityAdvectionType(Dropdown dropdown)
    {
        if (dropdown.value == 0) GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityAdvectionType = solver.VelocityAdvection.constant;
        else if (dropdown.value == 1) GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityAdvectionType = solver.VelocityAdvection.variable;
    }

    public void setDensityAdvection(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().densityAdvection = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().densityAdvection.ToString();
    }

    public void setVelocityAdvection(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityAdvection = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityAdvection.ToString();
    }

    public void setDensityDissipation(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().densityDissipation = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().densityDissipation.ToString();
    }

    public void setVelocityDissipation(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityDissipation = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().velocityDissipation.ToString();
    }

    public void conserveMass(Toggle toggle)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().conserveMass = toggle.isOn;
    }

    public void setSolverQualityFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().solverQuality = (int)slider.value;
        UpdateField(slider);
    }

    public void setSolverQualityFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().solverQuality = (int)Mathf.Clamp(Convert.ToInt32(field.text), solverQuality.x, solverQuality.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().solverQuality.ToString();
        UpdateSlider(field);
    }

    public void displayVelocity(Toggle toggle)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().displayVelocity = toggle.isOn;
    }

    public void setVelocityMultiplier(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().displayVelocityMultiplier = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().displayVelocityMultiplier.ToString();
    }

    public void setFlowDisplay(Dropdown dropdown)
    {
        if (dropdown.value == 0) GameObject.Find("wing_poly_geo").GetComponent<solver>().displayG = solver.DisplayG.none;
        else if (dropdown.value == 1) GameObject.Find("wing_poly_geo").GetComponent<solver>().displayG = solver.DisplayG.density;
        else if (dropdown.value == 2) GameObject.Find("wing_poly_geo").GetComponent<solver>().displayG = solver.DisplayG.pressure;
    }

    public void setFlowDisplayMultiplier(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().displayGridMultiplier = Mathf.Clamp(Convert.ToSingle(field.text), 0, 9999);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().displayGridMultiplier.ToString();
    }

    public void setAirfoilType(Dropdown dropdown)
    {
        if (dropdown.value == 0)
        {
            GameObject.Find("wing_poly_geo").GetComponent<solver>().airfoilType = solver.DisplayGrid.joukowsky;
            GameObject.Find("displayJoukowsky").transform.localScale = Vector3.one;
            GameObject.Find("xDisplacement").transform.localScale = Vector3.one;
            GameObject.Find("yDisplacement").transform.localScale = Vector3.one;
            GameObject.Find("maxCamber4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications4m").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5m").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications5m").transform.localScale = Vector3.zero;
            GameObject.Find("joukowsky").transform.localScale = Vector3.one;
            GameObject.Find("nacaFourDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitModifiedAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitModifiedAirfoil").transform.localScale = Vector3.zero;
        }
        else if (dropdown.value == 1)
        {
            GameObject.Find("wing_poly_geo").GetComponent<solver>().airfoilType = solver.DisplayGrid.nacaFourDigit;
            GameObject.Find("displayJoukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("xDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("yDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4").transform.localScale = Vector3.one;
            GameObject.Find("maxCamberPosition4").transform.localScale = Vector3.one;
            GameObject.Find("maxThickness4").transform.localScale = Vector3.one;
            GameObject.Find("maxCamber4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications4m").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5m").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications5m").transform.localScale = Vector3.zero;
            GameObject.Find("joukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitAirfoil").transform.localScale = Vector3.one;
            GameObject.Find("nacaFourDigitModifiedAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitModifiedAirfoil").transform.localScale = Vector3.zero;
        }
        else if (dropdown.value == 2)
        {
            GameObject.Find("wing_poly_geo").GetComponent<solver>().airfoilType = solver.DisplayGrid.nacaFourDigitModified;
            GameObject.Find("displayJoukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("xDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("yDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4m").transform.localScale = Vector3.one;
            GameObject.Find("maxCamberPosition4m").transform.localScale = Vector3.one;
            GameObject.Find("maxThickness4m").transform.localScale = Vector3.one;
            GameObject.Find("modifications4m").transform.localScale = Vector3.one;
            GameObject.Find("designLiftCoefficient5").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5m").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications5m").transform.localScale = Vector3.zero;
            GameObject.Find("joukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitModifiedAirfoil").transform.localScale = Vector3.one;
            GameObject.Find("nacaFiveDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitModifiedAirfoil").transform.localScale = Vector3.zero;
        }
        else if (dropdown.value == 3)
        {
            GameObject.Find("wing_poly_geo").GetComponent<solver>().airfoilType = solver.DisplayGrid.nacaFiveDigit;
            GameObject.Find("displayJoukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("xDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("yDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications4m").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5").transform.localScale = Vector3.one;
            GameObject.Find("maxCamberPosition5").transform.localScale = Vector3.one;
            GameObject.Find("camberType5").transform.localScale = Vector3.one;
            GameObject.Find("maxThickness5").transform.localScale = Vector3.one;
            GameObject.Find("designLiftCoefficient5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5m").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications5m").transform.localScale = Vector3.zero;
            GameObject.Find("joukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitModifiedAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitAirfoil").transform.localScale = Vector3.one;
            GameObject.Find("nacaFiveDigitModifiedAirfoil").transform.localScale = Vector3.zero;
        }
        else if (dropdown.value == 4)
        {
            GameObject.Find("wing_poly_geo").GetComponent<solver>().airfoilType = solver.DisplayGrid.nacaFiveDigitModified;
            GameObject.Find("displayJoukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("xDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("yDisplacement").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamber4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition4m").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness4m").transform.localScale = Vector3.zero;
            GameObject.Find("modifications4m").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5").transform.localScale = Vector3.zero;
            GameObject.Find("maxCamberPosition5").transform.localScale = Vector3.zero;
            GameObject.Find("camberType5").transform.localScale = Vector3.zero;
            GameObject.Find("maxThickness5").transform.localScale = Vector3.zero;
            GameObject.Find("designLiftCoefficient5m").transform.localScale = Vector3.one;
            GameObject.Find("maxCamberPosition5m").transform.localScale = Vector3.one;
            GameObject.Find("camberType5m").transform.localScale = Vector3.one;
            GameObject.Find("maxThickness5m").transform.localScale = Vector3.one;
            GameObject.Find("modifications5m").transform.localScale = Vector3.one;
            GameObject.Find("joukowsky").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFourDigitModifiedAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitAirfoil").transform.localScale = Vector3.zero;
            GameObject.Find("nacaFiveDigitModifiedAirfoil").transform.localScale = Vector3.one;
        }
    }

    public void setAngleOfAttackFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().angleOfAttack = slider.value;
        UpdateField(slider);
    }

    public void setAngleOfAttackFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().angleOfAttack = Mathf.Clamp(Convert.ToInt32(field.text), angleOfAttack.x, angleOfAttack.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().angleOfAttack.ToString();
        UpdateSlider(field);
    }

    public void setJoukowskyDisplay(Dropdown dropdown)
    {
        if (dropdown.value == 0) GameObject.Find("wing_poly_geo").GetComponent<solver>().displayJoukowsky = solver.DisplayJoukowsky.analyticVelocity;
        else if (dropdown.value == 1) GameObject.Find("wing_poly_geo").GetComponent<solver>().displayJoukowsky = solver.DisplayJoukowsky.analyticPressure;
        else if (dropdown.value == 2) GameObject.Find("wing_poly_geo").GetComponent<solver>().displayJoukowsky = solver.DisplayJoukowsky.numericSolver;
    }

    public void setXDisplacementFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().xDisplacement = slider.value;
        UpdateField(slider);
    }

    public void setXDisplacementFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().xDisplacement = Mathf.Clamp(Convert.ToSingle(field.text), xDisplacement.x, xDisplacement.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().xDisplacement.ToString();
        UpdateSlider(field);
    }

    public void setYDisplacementFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().yDisplacement = slider.value;
        UpdateField(slider);
    }

    public void setYDisplacementFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().yDisplacement = Mathf.Clamp(Convert.ToSingle(field.text), yDisplacement.x, yDisplacement.y);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().yDisplacement.ToString();
        UpdateSlider(field);
    }

    public void setMaxCamber4FromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4 = slider.value;
        UpdateField(slider);
    }

    public void setMaxCamber4FromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4 = Mathf.Clamp(Convert.ToSingle(field.text), 0, 25);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4.ToString();
        UpdateSlider(field);
    }

    public void setMaxCamberPosition4FromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4 = slider.value;
        UpdateField(slider);
    }

    public void setMaxCamberPosition4FromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4 = Mathf.Clamp(Convert.ToSingle(field.text), 0.5f, 9.5f);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4.ToString();
        UpdateSlider(field);
    }

    public void setMaxThickness4FromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4 = slider.value;
        UpdateField(slider);
    }

    public void setMaxThickness4FromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4 = Mathf.Clamp(Convert.ToSingle(field.text), 0, 30);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4.ToString();
        UpdateSlider(field);
    }

    public void setMaxCamber4mFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4m = slider.value;
        UpdateField(slider);
    }

    public void setMaxCamber4mFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4m = Mathf.Clamp(Convert.ToSingle(field.text), 0, 15);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamber4m.ToString();
        UpdateSlider(field);
    }

    public void setMaxCamberPosition4mFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4m = slider.value;
        UpdateField(slider);
    }

    public void setMaxCamberPosition4mFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4m = Mathf.Clamp(Convert.ToSingle(field.text), 0.5f, 9.5f);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition4m.ToString();
        UpdateSlider(field);
    }

    public void setMaxThickness4mFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4m = slider.value;
        UpdateField(slider);
    }

    public void setMaxThickness4mFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4m = Mathf.Clamp(Convert.ToSingle(field.text), 0, 30);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness4m.ToString();
        UpdateSlider(field);
    }

    public void setModifications4m(Dropdown dropdown)
    {
        if (dropdown.value == 0) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.Three;
        else if (dropdown.value == 1) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.Five;
        else if (dropdown.value == 2) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.ThirtyThree;
        else if (dropdown.value == 3) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.ThirtyFour;
        else if (dropdown.value == 4) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.ThirtyFive;
        else if (dropdown.value == 5) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.SixtyTwo;
        else if (dropdown.value == 6) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.SixtyThree;
        else if (dropdown.value == 7) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.SixtyFour;
        else if (dropdown.value == 8) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.SixtyFive;
        else if (dropdown.value == 9) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications4m = solver.Display4m.NinetyThree;
    }

    public void setDesignLiftCoefficient5FromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5 = slider.value;
        UpdateField(slider);
    }

    public void setDesignLiftCoefficient5FromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5 = Mathf.Clamp(Convert.ToSingle(field.text), 0.05f, 6.67f);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5.ToString();
        UpdateSlider(field);
    }

    public void setMaxCamberPosition5FromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5 = (int)slider.value;
        UpdateField(slider);
    }

    public void setMaxCamberPosition5FromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5 = (int)Mathf.Clamp(Convert.ToSingle(field.text), 1, 5);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5.ToString();
        UpdateSlider(field);
    }

    public void setCamberType5FromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5 = (int)slider.value;
        UpdateField(slider);
    }

    public void setCamberType5FromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5 = (int)Mathf.Clamp(Convert.ToSingle(field.text), 0, 1);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5.ToString();
        UpdateSlider(field);
    }

    public void setMaxThickness5FromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5 = slider.value;
        UpdateField(slider);
    }

    public void setMaxThickness5FromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5 = Mathf.Clamp(Convert.ToSingle(field.text), 0, 30);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5.ToString();
        UpdateSlider(field);
    }

    public void setDesignLiftCoefficient5mFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5m = slider.value;
        UpdateField(slider);
    }

    public void setDesignLiftCoefficient5mFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5m = Mathf.Clamp(Convert.ToSingle(field.text), 0.05f, 6.67f);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().designLiftCoefficient5m.ToString();
        UpdateSlider(field);
    }

    public void setMaxCamberPosition5mFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5m = (int)slider.value;
        UpdateField(slider);
    }

    public void setMaxCamberPosition5mFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5m = (int)Mathf.Clamp(Convert.ToSingle(field.text), 1, 5);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxCamberPosition5m.ToString();
        UpdateSlider(field);
    }

    public void setCamberType5mFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5m = (int)slider.value;
        UpdateField(slider);
    }

    public void setCamberType5mFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5m = (int)Mathf.Clamp(Convert.ToSingle(field.text), 0, 1);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().camberType5m.ToString();
        UpdateSlider(field);
    }

    public void setMaxThickness5mFromSlider(Slider slider)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5m = slider.value;
        UpdateField(slider);
    }

    public void setMaxThickness5mFromField(InputField field)
    {
        GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5m = Mathf.Clamp(Convert.ToSingle(field.text), 0, 30);
        field.text = GameObject.Find("wing_poly_geo").GetComponent<solver>().maxThickness5m.ToString();
        UpdateSlider(field);
    }

    public void setModifications5m(Dropdown dropdown)
    {
        if (dropdown.value == 0) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.Three;
        else if (dropdown.value == 1) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.Five;
        else if (dropdown.value == 2) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.ThirtyThree;
        else if (dropdown.value == 3) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.ThirtyFour;
        else if (dropdown.value == 4) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.ThirtyFive;
        else if (dropdown.value == 5) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.SixtyTwo;
        else if (dropdown.value == 6) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.SixtyThree;
        else if (dropdown.value == 7) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.SixtyFour;
        else if (dropdown.value == 8) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.SixtyFive;
        else if (dropdown.value == 9) GameObject.Find("wing_poly_geo").GetComponent<solver>().modifications5m = solver.Display5m.NinetyThree;
    }

    public void setEmissionType(Dropdown dropdown)
    {
        if (dropdown.value == 0)
        {
            GameObject.Find("inletEmission").transform.localScale = Vector3.one;
            GameObject.Find("mouseEmission").transform.localScale = Vector3.zero;
        }

        else if (dropdown.value == 1)
        {
            GameObject.Find("inletEmission").transform.localScale = Vector3.zero;
            GameObject.Find("mouseEmission").transform.localScale = Vector3.one;
        }
    }
}
