using UnityEngine;
using System;
using System.Linq;
using System.IO;
using System.Collections;
using System.Collections.Generic;

//  Código de inicialización de parámetros y variables y asignación de valores.

public static class Fresnel
{
    // Aquí se incluye todo lo previo al bucle de simulacion. Es el Setup del proyecto.

    /////////////////////////////////////
    //  Variables para el Setup
    /////////////////////////////////////
    public static int numeroperiodosint;
    public static double diajuliano1, angulodiario1, Ts1;
    public static double tiempomax, Ntotal, tactualizacion;
    public static double[] caudal1, Tentrada1, Tsalida1, I1, Tambiente1, minuto1, hora1;
    public static double[] Tf, Tm;

    /////////////////////////////////////
    //  Variables para Inicializacion
    /////////////////////////////////////
    public static double Latitud, oriplanta, Longitud, Eficienciamedia;
    public static int incrementodetiempo, tacteficiencias;
    public static double tint, diametrointeriortubo, diametroexteriortubo,
                        diametrointeriortubo2, diametroexteriortubo2;
    public static float Absortancia, Transmisividad, Reflectividad, Reflectividadsecundaria,
                        factorensuciamiento1, factorensuciamiento2;
    public static long tiemposimulacion;
    public static double Af, Am, G, L, pm, Cm, Am2, Af2, L2;
    public static int anyo, mes, dia, i;
    public static int numeroPartesDiscretasSistema = 254;
    public static int numeroPartesDiscretasTuboCaptador = 64;

    public static void Setup()
    {
        //  Primero inicializamos los datos del sistema del captador.
        Inicializacion();

        //  Ahora preparamos el sistema para operación.


        //--1--// Cargamos los datos experimentales

        ObtenerDatos();
        //  Aprovechamos que son variables públicas, no tenemos ni que recogerlas.


        //--2--// Interpolamos los datos de las medidas, para tener una base de muestras mayor

        numeroperiodosint = (int)Math.Ceiling(incrementodetiempo / tint);
        hora1       = Interpola(numeroperiodosint, hora1);
        caudal1     = Interpola(numeroperiodosint, caudal1);
        Tentrada1   = Interpola(numeroperiodosint, Tentrada1);
        Tsalida1    = Interpola(numeroperiodosint, Tsalida1);
        I1          = Interpola(numeroperiodosint, I1);
        Tambiente1  = Interpola(numeroperiodosint, Tambiente1);


        //--3--// Inicialización de las variables.

        diajuliano1 = FresnelSupport.CalculoDiaJuliano(anyo, mes, dia);

        angulodiario1 = FresnelSupport.CalculoAnguloDiario(diajuliano1);

        minuto1 = new double[hora1.Length];
        for(i = 0; i < hora1.Length; i++)
        {
            minuto1[i] = Math.Floor(Math.Abs(hora1[i] - Math.Truncate(hora1[i]))*60);
        }

        // Aplicamos, para inicializar el parámetro Ts1, la función, al primer elemento de los datos experimentales
        Ts1 = FresnelSupport.CalculoHoraSolar(hora1[0], minuto1[0], mes, angulodiario1);

        Tf = new double[numeroPartesDiscretasSistema];   // Tiene la longitud de las partes de fluido que consideramos en el sistema.
        for(i = 0; i < numeroPartesDiscretasSistema; i++)
        {
            Tf[i] = Tentrada1[0];   // Inicializar la temperatura del fluido
        }

        Tm = new double[Tf.Length];
        for(i = 0; i < numeroPartesDiscretasSistema; i++)
        {
            Tm[i] = Tf[i] + 7;    // Inicializar la temperatura del metal
        }


        //--4--// Cálculos necesarios para el funcionamiento del programa.

        //  Cálculo del número de iteraciones, correspondientes a los tiempos de
        //  simulación e integración que se han definido en el fichero de
        //  configuración. También se calcula el tiempo máximo de simulación.
        //-> tiempomax es el número máximo de iteraciones que hará el bucle principal.
        // Es así porque se corresponde con la longitud del archivo de datos (podríamos haber
        // puesto cualquier otra de las variables que contiene, aquí usamos "Tentrada1").
        tiempomax = Tentrada1.Length;
        
        //-> Ntotal es el número de iteraciones del bucle principal de simulación,
        // calculado en base a los tiempos que hemos indicado en Inicializacion.
        // Si es mayor que tiempomax, se sobreescribe, pues no tendríamos datos para hacer
        // la simulación.
        Ntotal = Math.Floor(Fresnel.tiemposimulacion / Fresnel.tint);
        if(Ntotal > tiempomax) Ntotal = tiempomax;

        // Cada cuántas muestras se actualizan las eficiencias
        //-> tactualizacion marca la frecuencia de actualización de FactorSombra, que aglutina bastantes cálculos
        // sobre la eficiencia geométrica y las sombras posibles en el captador.
        tactualizacion = Math.Ceiling(Fresnel.tacteficiencias / Fresnel.tint);

        Debug.Log("Fresnel se ha cargado");
    }

    private static void Inicializacion()
    {
        Debug.Log("Los datos del sistema se han cargado");
        
        // Introducir fecha correspondiente a los datos que se usarán para la simulación
        // En realidad estos datos tenemos que cargarlos de un fichero, estas declaraciones no me sirven

        // anyo  = 2009;
        // mes   = 06;
        // dia   = 16;
        
        anyo  = 2009;
        mes   = 06;
        dia   = 29;

        // anyo  = 2009;
        // mes   = 07;
        // dia   = 06;

        // anyo  = 2009;
        // mes   = 07;
        // dia   = 10;

        // anyo  = 2009;
        // mes   = 08;
        // dia   = 17;

        // anyo  = 2009;
        // mes   = 08;
        // dia   = 22;
        


        // ------------------->>PARTE CONFIGURABLE DEL FICHERO>>--------------------
        //__________________________________________________________________________
        // Tiempo de integración, Tiempo de muestreo de sensores, Simulación,
        // ---------------------Representación en segundos--------------------------
        // Tiempo de integración
        tint = 0.25;

        // Tiempo de muestreo de sensores. 
        // Si se dan cada 15 minutos son 900 segundos.
        incrementodetiempo = 60; 

        // Tiempo de simulación en segundos
        tiemposimulacion = 7*3600; 

        // El tiempo para actualizar las eficiencias en segundos.
        tacteficiencias = 10;


        // Asignacion de valor de variables dimensionales

        // Tuberia 1 [m]
        diametrointeriortubo = 0.068; //[m]
        diametroexteriortubo = 0.071; //[m]

        // Tuberia 2 [m]
        diametrointeriortubo2 = 0.076986659489532;    //[m]
        diametroexteriortubo2 = 0.16;                 //[m]

        // Apertura del colector [m]
        G = 0.5*11; 

        // Características del metal DIN 1.4404
        pm = 8027;    //Densidad del metal (DIN 1.4404) [kg/m3]
        Cm = 500;     //Calor específico del metal      [J/kg°C] 

        // MODELO ÓPTICO DEL CAPTADOR SOLAR FRESNEL
        //__________________________________________________________________________
        //  Variables de orientación y posición:Latitud y longitud 
        //  y orientacion e inclinación del campo

        Longitud    = -6;
        Latitud     = 37.41;
        oriplanta   = (12 + 3/600 + 11/3600)*Math.PI/180;
            // ¿Qué sentido tiene esta transformación?
            // Entiendo que se debe a que la planta tiene desviacion de 
            // 12º 3' 1'' de la direccion sur (pero entonces no entiendo el '11' en la formula, ni el 600 dividiendo)

        // Dimensiones y propiedades del tubo receptor.
        //__________________________________________________________________________
        //   Definición de reflectividad y factores de eficiencia.
        //   Se ha supuesto que si está limpio factorensuciamiento1=0.9 
        //   y factorensuciamiento2=0.9, y si está sucio 0.8 para ambos.
        Absortancia     = 0.94f;
        Transmisividad  = 0.96f;
        Reflectividad   = 0.9f*0.9f;

        Reflectividadsecundaria = 0.77f;
        factorensuciamiento1    = 0.8f;
        factorensuciamiento2    = 0.8f;

        // -------->>AQUÍ COMIENZA LA PARTE NO CONFIGURABLE DEL FICHERO<<----------
        //                          CÁLCULOS INICIALES.
        //__________________________________________________________________________

        // Cálculo de las eficiencias
        Eficienciamedia = Absortancia*Reflectividad*Transmisividad*Reflectividadsecundaria
                            *factorensuciamiento1*factorensuciamiento2;
        // Para Tubería 1
        // Diametro propio de la tuberia [m]        
        L = Math.PI*diametrointeriortubo; 

        // Áreas transversales del fluido[m2]
        Af = Math.PI*Math.Pow(diametrointeriortubo, 2)/4;
        Am = Math.PI*(Math.Pow(diametroexteriortubo, 2) - Math.Pow(diametrointeriortubo, 2))/4;

        // Para Tubería 2
        // Diametro propio de la tuberia 2[m]
        L2 = Math.PI*diametrointeriortubo2; 

        // Áreas transversales del fluido [m2]
        Af2 = Math.PI*Math.Pow(diametrointeriortubo2, 2)/4;
        Am2 = Math.PI*(Math.Pow(diametroexteriortubo2, 2) - Math.Pow(diametrointeriortubo2, 2))/4;//(m2)
    }

    private static void ObtenerDatos()
    {
        string[,] DatosCSV = CSVReader.LeerCSV(@"C:\Users\Pablo\OneDrive - UNIVERSIDAD DE SEVILLA\TFG\0_Integracion\Assets\DataFiles\datosMedidos.csv");

        bool success;
        int longitudBucle = DatosCSV.GetLength(0);

        hora1      = new double[DatosCSV.GetLength(0)];
        caudal1    = new double[DatosCSV.GetLength(0)];
        Tentrada1  = new double[DatosCSV.GetLength(0)];
        Tsalida1   = new double[DatosCSV.GetLength(0)];
        I1         = new double[DatosCSV.GetLength(0)];
        Tambiente1 = new double[DatosCSV.GetLength(0)];

        for(int contador = 0; contador < longitudBucle; contador++)
        {
            success = Double.TryParse(DatosCSV[contador, 0], out hora1[contador]);      // h
            if(!success) Debug.Log("Error al convertir hora a formato double, en el indice: " + contador);
            success = Double.TryParse(DatosCSV[contador, 1], out caudal1[contador]);    // m3/s
            if(!success) Debug.Log("Error al convertir caudal a formato double, en el indice: " + contador);
            success = Double.TryParse(DatosCSV[contador, 2], out Tentrada1[contador]);  // ºC
            if(!success) Debug.Log("Error al convertir Tª entrada a formato double, en el indice: " + contador);
            success = Double.TryParse(DatosCSV[contador, 3], out Tsalida1[contador]);   // ºC
            if(!success) Debug.Log("Error al convertir Tª salida a formato double, en el indice: " + contador);
            success = Double.TryParse(DatosCSV[contador, 5], out I1[contador]);         // W/m2
            if(!success) Debug.Log("Error al convertir irradiancia a formato double, en el indice: " + contador);
            success = Double.TryParse(DatosCSV[contador, 6], out Tambiente1[contador]); // ºC
            if(!success) Debug.Log("Error al convertir Tª ambiente a formato double, en el indice: " + contador);
        }
    }

    private static double[] Interpola(int numeroMuestras, double[] datosOrigen)
    {
        int tam = datosOrigen.Length;
        int tamvect = numeroMuestras;
        int n = tamvect;

        // Primera iteración fuera del bucle para acumular
        double[] vectorsalida = new double[tam * tamvect];

        for(int i = 1; i <= tamvect; i++)
        {
            vectorsalida[i - 1] = datosOrigen[0] + ((datosOrigen[1] - datosOrigen[0]) / numeroMuestras) * i;
        }

        for(int i = 1; i < tam; i++)
        {
            for(int j = 1; j <= tamvect; j++)
            {
                if(i >= tam - 1)
                {
                    vectorsalida[n] = datosOrigen[i - 1] + ((datosOrigen[i] - datosOrigen[i - 1]) / numeroMuestras) * j;
                }
                else vectorsalida[n] = datosOrigen[i] + ((datosOrigen[i + 1] - datosOrigen[i]) / numeroMuestras) * j;
                n++;
            }
        }

        return vectorsalida;
    }

}