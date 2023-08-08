using UnityEngine;
using System;
using System.Linq;
using System.IO;
using System.Collections;
using System.Collections.Generic;

//  Código de inicialización de parámetros y variables y asignación de valores.

public static class Fresnel
{
    // Aqui tengo que incluir todo lo previo al bucle de simulacion. Va a ser el Setup del proyecto.

    /////////////////////////////////////
    //  Variables para el Setup
    /////////////////////////////////////
    public static int numeroperiodosint, k1;
    public static double    diajuliano1, angulodiario1, factorsombra1, dia1, mes1, anyo1,
                            irradiancia1, tambiente1, q1, Ts1;
    public static double tiempomax, Ntotal, tactualizacion, incrtiempo;
    public static double[] caudal1, Tentrada1, Tsalida1, I1, Tambiente1, minuto1, hora1;
    public static double[] Tf, Tm;
    public static string[,] DatosExperimentales;

    /////////////////////////////////////
    //  Variables para Inicializacion
    /////////////////////////////////////
    public static double Latitud, oriplanta, Longitud, Eficienciamedia, horacomienzo;
    public static int incrementodetiempo, tacteficiencias, tiempoderepresentacion;
    public static double tint, dt, reext, diametrointeriortubo, diametroexteriortubo,
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

        mes1 = mes; // Esto mejor lo quitamos, directamente eliminamos las variables: anyo1, mes1, dia1;
                    // puesto que ya contamos con las simples, porque hemos incluido "DatosSistema" aqui.
        diajuliano1 = FresnelSupport.CalculoDiaJuliano(anyo, mes, dia);

        angulodiario1 = FresnelSupport.CalculoAnguloDiario(diajuliano1);

        minuto1 = new double[hora1.Length];
        for(i = 0; i < hora1.Length; i++)
        {
            // Rematar este cálculo en busca de los dos decimales, para que tenga sentido? // ***
            minuto1[i] = Math.Floor(Math.Abs(hora1[i] - Math.Truncate(hora1[i]))*60);
        }

        // Aplicamos, para inicializar el parametro Ts1, la funcion, al primer elemento de los datos experimentales
        Ts1 = FresnelSupport.CalculoHoraSolar(hora1[0], minuto1[0], mes, angulodiario1);

        factorsombra1 = FresnelSupport.eficienciaGeoYSombras(angulodiario1, Ts1);

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
        // puesto cualquier otra de las variables que contiene).
        tiempomax = Tentrada1.Length;
        
        //-> Ntotal es el número de iteraciones del bucle principal de simulación,
        // calculado en base a los tiempos que hemos indicado en Inicializacion.
        // Si es mayor que tiempomax, se sobreescribe, pues no tendríamos datos para hacer
        // la simulación. *** Revisar esto por si tiene sentido anularlo y ceñirse al número de muestras que tenemos.
        Ntotal = Math.Floor(Fresnel.tiemposimulacion / Fresnel.tint);
        if(Ntotal > tiempomax) Ntotal = tiempomax;

        // Cada cuantas muestras se actualizan las eficiencias
        //-> tactualizacion marca la frecuencia de actualización de FactorSombra, que aglutina bastantes cálculos
        // sobre la eficiencia geométrica y las sombras posibles en el captador.
        tactualizacion = Math.Ceiling(Fresnel.tacteficiencias / Fresnel.tint); // *** Revisarlo y corregir

        // Cada cuantas muestras hay que leer los datos del fichero.
        //-> Ahora mismo no lo usamos para nada. A expensas de aclarar cuál era el objetivo.
        incrtiempo = Math.Ceiling(Fresnel.incrementodetiempo / Fresnel.tint); // *** Revisarlo y corregir

        // Para representacion grafica
        // horarep = 11;   // *** eliminarlo si no sirve

        k1 = 0; // k es un parámetro que determinamos en SetUp, por si quisiéramos empezar en una muestra determinada.
                // De todas formas, para realizar eso con éxito, deberíamos parametrizarlo de manera que el tiempo fuese acorde.

        Debug.Log("Fresnel se ha cargado");
    }

    private static void Inicializacion()
    {
        Debug.Log("Los datos del sistema se han cargado");
        
        // Introducir fecha correspondiente a los datos que se usaran para la simulacion
        // En realidad estos datos tenemos que cargarlos de un fichero, estas declaraciones no me sirven
        anyo  = 2009;   //  System.DateTime.Now.Year;
        mes   = 06;     //  System.DateTime.Now.Month;
        dia   = 29;     //  System.DateTime.Now.Day;
        // hora  = System.DateTime.Now.TimeOfDay.Hours;
        // minuto = System.DateTime.Now.TimeOfDay.Minutes;

        // Tiempo de integración
        dt      = 0.5;
        reext   = 0.0225;

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
        tiempoderepresentacion = 10; //cada minuto

        // El tiempo para actualizar las eficiencias en segundos.
        tacteficiencias = 10; 

        // hora de comienzo en decimal
        // *** pendiente de corregir, creo que no queremos indicar una hora de comienzo de las medidas, simplemente
        //    tomaremos la del fichero
        horacomienzo = 11; // A ver esto de donde sale

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
            // ¿Que sentido tiene esta transformacion?
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

        // recoger las eficiencias de ficheros
        Eficienciamedia = Absortancia*Reflectividad*Transmisividad*Reflectividadsecundaria
                            *factorensuciamiento1*factorensuciamiento2;
        // Para Tuberia 1
        // Diametro propio de la tuberia [m]        
        L = Math.PI*diametrointeriortubo; 

        // Areas transversales del fluido[m2]
        Af = Math.PI*Math.Pow(diametrointeriortubo, 2)/4; 
        Am = Math.PI*(Math.Pow(diametroexteriortubo, 2) - Math.Pow(diametrointeriortubo, 2))/4;

        // Para Tuberia 2
        // Diametro propio de la tuberia 2[m]
        L2 = Math.PI*diametrointeriortubo2; 

        // Areas transversales del fluido [m2]
        Af2 = Math.PI*Math.Pow(diametrointeriortubo2, 2/4);
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
        int tam = datosOrigen.Length - 1;
        int tamvect = numeroMuestras - 1;
        int n = tamvect;

        // primera iteracion fuera del bucle para acumular
        double[] vectorsalida = new double[tam * tamvect];

        for(int i = 1; i <= tamvect; i++)
        {
            vectorsalida[i - 1] = datosOrigen[0] + ((datosOrigen[1] - datosOrigen[0]) / numeroMuestras) * i;
            // Debug.Log("VectorSalida: "+vectorsalida[i-1]);
        }

        for(int i = 1; i < tam; i++)
        {
            for(int j = 1; j <= tamvect; j++)
            {
                vectorsalida[n] = datosOrigen[i] + ((datosOrigen[i + 1] - datosOrigen[i]) / numeroMuestras) * j;
                // Debug.Log("VectorSalida: "+vectorsalida[n]);
                n++;
            }
        }

        return vectorsalida;
    }

}