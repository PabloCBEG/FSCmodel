using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;
// using CsvHelper;
// using CsvHelper.Configuration.Attributes;
// using Plotly.NET.CSharp;

// #r "nuget: Plotly.NET, 4.2.0";

public class Simulacion : MonoBehaviour
{
    double[,] TemperaturasSimulacion;
    double[] TemperaturaFluidoSimulacion, TemperaturaMetalSimulacion;
    // double[]        TemperaturaSalidaFluido, TemperaturaTuberia, TemperaturaEntradaFluido,
    //                 CaudalSalida, IrradianciaSalida;
    double   IrradianciaSimulacion, PorcentajeTuboSimulacion, TiempoSolarSimulacion, FactorSombraSimulacion,
                    TemperaturaAmbienteSimulacion, qSimulacion, EficienciaMediaSimulacion,
                    minutoSimulacion, horaSimulacion, RepeticionesSubcicloSimulacion;
    int j, i, n, Ntotal, contador10;

    string filePath = @"C:\Users\Pablo\OneDrive - UNIVERSIDAD DE SEVILLA\TFG\FSCmodel\FSCmodelUnity\Assets\Scripts\Debugging_data_060923.txt";

    void Awake()
    {
        Debug.Log("Simulacion ha empezado");

        // Modifica el tiempo fijado de repetición del bucle. Para que tuviera sentido del todo,
        // habría que ejecutar los cálculos cada 10 iteraciones del bucle (así, el tiempo de
        // muestreo es realmente 0.25 s, como es en el modelo en matlab). De todas formas, los
        // resultados son aceptables sin llevar a cabo esta corrección sugerida.
        Time.fixedDeltaTime = 0.025f;

        
        // Se pone en marcha el sistema:
        Fresnel.Setup();

        TemperaturaFluidoSimulacion = Fresnel.Tf;
        TemperaturaMetalSimulacion  = Fresnel.Tm;

        // A hora y minuto de la simulación se le dan los valores del primer elemento de las muestras,
        // porque a partir del primer paso, el mismo bucle de simulación las calculará.
        minutoSimulacion            = Fresnel.minuto1[0];
        // Atencion al espíritu o la intención de la declaración en el modelo original: no es lo que estamos
        // haciendo: (original:) x = find(abs(hora-horacomienzo)==min(abs(hora-horacomienzo))); hora = hora(x);
        horaSimulacion              = Fresnel.hora1[0];

        EficienciaMediaSimulacion = Fresnel.Eficienciamedia;
        j = 0;

        Ntotal = (int)Fresnel.Ntotal;
        

        // TemperaturaSalidaFluido     = new double[Ntotal];
        // TemperaturaTuberia          = new double[Ntotal];
        // TemperaturaEntradaFluido    = new double[Ntotal];
        // CaudalSalida                = new double[Ntotal];
        // IrradianciaSalida           = new double[Ntotal];

        contador10 = 1;
        
    }

    void FixedUpdate()
    {
        TemperaturaFluidoSimulacion[0]  = Fresnel.Tentrada1[j]; //  Inicializamos el primer valor (la temperatura de entrada del fluido)
        IrradianciaSimulacion           = Fresnel.I1[j];        //  Es constante para todo el sistema
        TemperaturaAmbienteSimulacion   = Fresnel.Tambiente1[j];//  Es constante para todo el sistema
        qSimulacion                     = Fresnel.caudal1[j];   //  Es constante para toda la tubería

        PorcentajeTuboSimulacion = 1;
        // Si la temperatura de salida es mayor que 180ºC
        if(TemperaturaFluidoSimulacion[Fresnel.numeroPartesDiscretasTuboCaptador] > (double)180) PorcentajeTuboSimulacion = 0;

        if(qSimulacion < 0) qSimulacion = 0;

        if(j % Fresnel.tactualizacion == 0)
        {
            minutoSimulacion = minutoSimulacion + (double)10 / (double)60;
            if(minutoSimulacion >= (double)60)
            {
                horaSimulacion++;
                if(horaSimulacion == 24) horaSimulacion = 0;
                minutoSimulacion -= (double)60;
            }

            TiempoSolarSimulacion = FresnelSupport.CalculoHoraSolar(horaSimulacion, minutoSimulacion, Fresnel.mes, Fresnel.angulodiario1);

            FactorSombraSimulacion = FresnelSupport.eficienciaGeoYSombras(Fresnel.angulodiario1, TiempoSolarSimulacion);
        }

        TemperaturasSimulacion = FresnelSupport.CalculoTemperatura( TemperaturaMetalSimulacion,
                                                                    TemperaturaFluidoSimulacion,
                                                                    TemperaturaAmbienteSimulacion,
                                                                    IrradianciaSimulacion,
                                                                    FactorSombraSimulacion,
                                                                    qSimulacion,
                                                                    EficienciaMediaSimulacion,
                                                                    PorcentajeTuboSimulacion);

        for(int indice = 0; indice < TemperaturasSimulacion.Length/2; indice++)
        {
            // El vector conjunto tiene longitud 2xlongitud de uno de los vectores de temperaturas
            TemperaturaMetalSimulacion[indice] = TemperaturasSimulacion[0, indice];
            TemperaturaFluidoSimulacion[indice] = TemperaturasSimulacion[1, indice];
        }

        // Store data for later plotting
        // In each iteration, append the value to the file
        File.AppendAllText(filePath, TemperaturaFluidoSimulacion[64].ToString() + Environment.NewLine);

        j++;

        if(j > Ntotal) UnityEditor.EditorApplication.isPlaying = false;
    }
}