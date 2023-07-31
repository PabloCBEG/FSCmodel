using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class SolarMovementSupport : MonoBehaviour
{
    // Parameters
    const int numberOfMirrorsInLine = 16;
    const int fresnelCollectorHeight = 4;   //metros
    public const float variacionMovimientoSolar = 0.5f;  //grados
    //public static double anguloEspejo;
    static int flagSol = 0;

    // Struct to save the position of the sun and the mirrors
    public struct FresnelSystem
    {
        public Transform Sun;
        public Transform[,] Mirrors;

    }

    // Function that find every mirror and the position of the sun
    public static FresnelSystem FindFresnelObjects()
    {
        GameObject[] aux;
        int i = 0, j = 0;

        FresnelSystem fresnelSystem;
        fresnelSystem.Sun = GameObject.Find("Sun").transform;
        aux = GameObject.FindGameObjectsWithTag("Mirror");
        fresnelSystem.Mirrors = new Transform[aux.Length/numberOfMirrorsInLine, numberOfMirrorsInLine];
        
        foreach(GameObject mirror in aux)
        {
            fresnelSystem.Mirrors[i, j] = mirror.transform;
            j++;
            if(j == numberOfMirrorsInLine)
            {
                i++;
                j = 0;
            }
        }

        return fresnelSystem;
    }

    public static double MirrorAngle(int rowNumber, float anguloRotacionAnterior)
    {
        // La funcion recibe el numero de fila del espejo, sabiendo que empezamos a contar por el espejo que en Unity
        // se corresponde con Fila1
        double anguloEspejo, h, angulo1;

        int coeficiente = Math.Abs(6 - rowNumber);
        
        // Calculo de la longitud del vector que une el centro del espejo con el tubo colector
        h = Math.Sqrt(Math.Pow(fresnelCollectorHeight, 2) + Math.Pow(coeficiente*0.7, 2));

        if(Math.Asin(coeficiente*0.7 / h) <= 90)
            angulo1 = 180*Math.Asin(coeficiente*0.7 / h)/Math.PI;
        else
            angulo1 = 180 - 180*Math.Asin(coeficiente*0.7 / h)/Math.PI;

        double angulo2 = 180 - 90 - angulo1;
        double angulo3 = 90 + angulo1;

        // Antes de hallar el angulo que debe girar el espejo, necesito saber qué angulo forma el sol
        double anguloSol = GameObject.Find("Sun").transform.localRotation.eulerAngles.x;
        
        if(flagSol == 1)    if(anguloSol < 90)  anguloSol = 180 - anguloSol;
        
        if(rowNumber > 5)   anguloEspejo = (angulo2 + anguloSol)/2;
        else                anguloEspejo = (angulo3 + anguloSol)/2;
        
        if(anguloSol == 90) flagSol = 1;

        return anguloEspejo;
    }

}

