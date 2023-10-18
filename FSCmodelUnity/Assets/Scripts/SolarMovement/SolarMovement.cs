using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class SolarMovement : MonoBehaviour
{
    // Struct to save the position of the sun and the mirrors
    SolarMovementSupport.FresnelSystem fresnelSystem;
    float rotacionActualSol, rotacionActualEspejo, anguloRotacionEspejo, anguloRotacionAnterior;
    float rotacionYplanta, rotacionZplanta;

    // Awake to find every mirror and the position of the sun
    void Awake()
    {
        fresnelSystem = SolarMovementSupport.FindFresnelObjects();
        rotacionYplanta = (float)GameObject.Find("Planta").transform.localRotation.eulerAngles.y;
        rotacionZplanta = (float)GameObject.Find("Planta").transform.localRotation.eulerAngles.z;
    }

    // Fixed Update to sims every instant the movement of the sun and the automatic control of every mirror
    void FixedUpdate()
    {
        fresnelSystem.Sun.Rotate(SolarMovementSupport.variacionMovimientoSolar, 0, 0);
        
        rotacionActualSol = (float)GameObject.Find("Sun").transform.localRotation.eulerAngles.x;

        // Lo hacemos por espejos, aunque para lo que estamos haciendo daría lo mismo hacerlo por filas,
        // con idea de implementarlo con más precisión en un futuro.
        for(int i = 0; i < fresnelSystem.Mirrors.GetLength(0); i++)
        {
            anguloRotacionEspejo = (float)SolarMovementSupport.MirrorAngle(i + 1, anguloRotacionAnterior);
            rotacionActualEspejo = (float)fresnelSystem.Mirrors[i, 0].transform.localRotation.eulerAngles.x;

            // Recorremos todas las filas y todas las columnas para girar los espejos
            for(int j = 0; j < fresnelSystem.Mirrors.GetLength(1); j++)
                fresnelSystem.Mirrors[i, j].transform.rotation = Quaternion.Euler(anguloRotacionEspejo, rotacionYplanta, rotacionZplanta);

            // Debug.Log("RotacionEspejo: "+anguloRotacionEspejo);

            anguloRotacionAnterior = anguloRotacionEspejo;
        }

        if(GameObject.Find("Sun").transform.localRotation.eulerAngles.x > 180) UnityEditor.EditorApplication.isPlaying = false;
    }
}