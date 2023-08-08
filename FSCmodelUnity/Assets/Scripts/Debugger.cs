
using UnityEngine;
using System;
using System.Linq;
using System.IO;
using System.Collections;
using System.Collections.Generic;

public class Debugger : MonoBehaviour
{
    private static string[,] DatosCSV;
    private static double[] Factorsombra1, Porcentajetubo1, TactualFluido1, TactualMetal1;
    static int counter;
    static double[,] TemperaturasSimulacion;
    private static string[] FactorsombraCSV, PorcentajetuboCSV, TactualFluidoCSV, TactualMetalCSV;

    void Start()
    {
        Fresnel.Setup();
        ObtenerDatos();
        counter = 0;
    }

    void FixedUpdate()
    {
        TemperaturasSimulacion = CalculoTemperatura( TactualMetal1[counter],
                                                                    TactualFluido1[counter],
                                                                    Fresnel.Tambiente1[counter],
                                                                    Fresnel.I1[counter],
                                                                    Factorsombra1[counter],
                                                                    Fresnel.caudal1[counter],
                                                                    Fresnel.Eficienciamedia,
                                                                    Porcentajetubo1[counter]);

        counter++;
    }

    public static double[,] CalculoCaracteristicasDelFluido(double[] Taceite, double Caudal1)
    {
        // Cálculo de las características físicas y termodinámicas del fluido, así *
        // como el coeficiente de calor por convección metal-fluido de cada trozo  *
        // en que se ha dividido el lazo de colectores.                            *
        // 
        // El cálculo se hace mediante ecuaciones ajustadas por mínimos cuadrados  *
        // con datos del fabricante.
        //                                                                         *
        // Entradas a la función:                                                  *
        //    - Temperatura del fluido [ºC]
        //    - caudal                 [m3/s]
        //                                                                         *
        // Salidas de la función:                                                  *
        //       - densidad pf [kg/m3]                                             *
        //       - Calor específico Cf [J/(kgºC)] 
        //       - Conductividad térmica Kf [W/(mºC)]
        //       - viscosidad mu [mPa*sec]
        //       - número de Prandtl Hv [adimensional]
        //       - Coeficiente de transmisión de calor por convección 
        //         metal-fluido Ht [W/(m2ºC)]

        //Caudal es un escalar, pero para cada muestra en cada momento, tiene un valor
        //(Tiene sentido, porque el caudal es constante a lo largo de la tuberia)
        //[ReadOnly] public double Caudal1;
        int tamanyoT;
        int l;
        double[] T, pf, Cf, Hv, Ht;
        double[,] resultado;

        //*** Por ahora tomamos la primera muestra de Caudal para los calculos;
        //*** mas adelante desarrollaremos para iterar todo.
        // Taceite = GetComponent<CalcTemp>().TactualFluido_CalcTemp;
        // Caudal1 = GetComponent<CalcTemp>().Caudal_CalcTemp;
        T = Taceite; // Cuidado porque en C# los arrays no se si se copian de un golpe o tengo que ir componente a componente
        tamanyoT = T.Length;

        // Debug.Log("Longitud de TemperaturaFluido en CalcTemp:"+tamanyoT);
        pf = new double[tamanyoT];
        Cf = new double[tamanyoT];
        Hv = new double[tamanyoT];
        Ht = new double[tamanyoT];
        resultado = new double[3, tamanyoT];

        // Aplicación de las ecuaciones
        for(l = 0; l < tamanyoT; l++)
        {
            pf[l] = - 0.002549807548*Math.Pow(T[l], 2)
                    - 0.202618731192*T[l]
                    + 1003.917572929273;
            Cf[l] = 0.000000516739*Math.Pow(T[l], 4)
                    - 0.000156862606*Math.Pow(T[l], 3)
                    + 0.027679455287*Math.Pow(T[l], 2)
                    - 1.626420964922*T[l]
                    + 4207.403959277281; 
            Hv[l] = 0.000000001338864e5*Math.Pow(T[l], 4)
                    - 0.000000778905228e5*Math.Pow(T[l], 3)
                    + 0.000187251179615e5*Math.Pow(T[l], 2)
                    - 0.025731157369768e5*T[l]
                    + 4.108383857418933e5;
            Ht[l] = Hv[l]*Math.Pow((Caudal1/3600), 0.8);

            resultado[0, l] = pf[l];
            resultado[1, l] = Cf[l];
            resultado[2, l] = Ht[l];
        }

        return resultado;
    }

    public static double[] CalculoPerdidasMetal(double[] T, double Tambiente)
    {
        // Cálculo del coeficiente de pérdidas metal-ambiente por m2 de superficie m
        // para cada trozo de tubo.
        //                                                                         *
        //  Entradas a la función:                                                 *
        //       - Temperatura metal [ºC]                                          *
        //       - Temperatura ambiente ºC]                                        *
        //                                                                         *
        // Salidas de la función:                                                  *
        //       - Coeficiente de pérdidas metal-ambiente [W/(m2ºC)]]              *
        //       - Pérdidas totales en todo el lazo       [KW]                     *

        double Sup; // Superficie metalica susceptible de perder energia en forma de calor
        int tamanyoT_PerdidasMetal, m;
        //private double[] T;  // Me parece una barbaridad declararlo en double,
                                                // pero no estoy seguro del espacio que necesita. Va a ocupar
                                                // muchisima memoria.
        double[] Hl;

        //*** Por ahora cogemos solo la primera muestra de Tambiente, luego desarrollaremos
        //*** el codigo para iterar todas.

        //Tambiente = GetComponent<CalcTemp>().Tambiente_CalcTemp;
        //tamanyoT_PerdidasMetal = GetComponent<CalcTemp>().TactualPerdidas.Length;
        //T = GetComponent<CalcTemp>().TactualPerdidas;
        tamanyoT_PerdidasMetal = T.Length;
        Hl = new double[tamanyoT_PerdidasMetal];

        Sup = 64*11*0.5;

        // Calculamos por metro cuadrado de espejo para perdidas de tubo             
        for(m = 0; m < tamanyoT_PerdidasMetal; m++)
        {
            Hl[m] = (4.5659247191149893E-1/Sup)*(T[m] - Tambiente) - 0.01062045206593238E+2/Sup;
            //Debug.Log("Hl: "+Hl[m]);
        }
        //Debug.Log("Hl: "+Hl[30]);

        //Hl(i)= 0.040251565842975*(T(i)-Tambiente)+ 3.528725247698581;
        // 0.018977269846620
        // 7.129191838532097;0.310522
        // Metros 0-64 de la tuberia se corresponden al captador;
        // Metros 65 en adelante, al resto de la tuberia, intercambiador de calor...
        // Debug.Log("Tamanyo Hl en perdidasMetal:"+Hl.Length);
        // Debug.Log("Tamanyo T en perdidasMetal:"+T.Length);
        for(m = 64; m <= 164; m++)
        {
            Hl[m] = 0.0140251565842975*(T[m] - Tambiente)*0.09 + 7.328725247698581*0.09;
        }
        for(m = 165; m<tamanyoT_PerdidasMetal; m++)
        {
            Hl[m] = 0*(T[m] - Tambiente) + 0.290223*5;
        }

        return Hl;
    }

    public static double[,] CalculoTemperatura( double[] TactualMetal_CalcTemp,
                                                double[] TactualFluido_CalcTemp,
                                                double Tambiente_CalcTemp,
                                                double I_CalcTemp,
                                                double Factorsombra_CalcTemp,
                                                double Caudal_CalcTemp,
                                                double Eficienciamedia_CalcTemp,
                                                double Porcentajetubo)
    // CalculoTemperatura recibe los parametros del sistema para cada instante de simulacion
    // Los dos primeros tienen un valor para cada discretizacion del tubo. Los demas, son escalares:
    // los suponemos constantes en toda la longitud del sistema.
    {
        // I en realidad es un array con valor para cada muestra.
        double Af, Af2, Am, Am2, G, L, L2, pm, Cm, tint_CalcTemp, q;
        double[] incx, Hl_local, pf_local, Cf_local, Ht_local, TfluidoSinCorregir;
        int n;
        double[] Tfluido_CalcTemp, Tmetal_CalcTemp, TactualPerdidas, vectorradiacion;
        double[,] aux2;
        double[,] TemperaturasCalcTemp;

        Af  = Fresnel.Af;  Af2 = Fresnel.Af2; Am  = Fresnel.Am;  Am2 = Fresnel.Am2;
        G   = Fresnel.G;   L = Fresnel.L;     L2 = Fresnel.L2;
        pm  = Fresnel.pm;  Cm  = Fresnel.Cm;
        tint_CalcTemp = Fresnel.tint;

        incx = new double[TactualFluido_CalcTemp.Length];
        for(n = 0; n < TactualFluido_CalcTemp.Length; n++)
        {
            if(n < 64) incx[n] = 1;
            else incx[n] = 0.9;
        }

        // Características del fluido  y coeficientes de pérdidas y transmisión de 
        // calor correspondientes a cada trozo de metal y fluido

        TactualPerdidas = new double[TactualFluido_CalcTemp.Length];
        for(n = 0; n < 64; n++)
        {
            TactualPerdidas[n] = TactualMetal_CalcTemp[n];
        }
        for(n = 64; n < TactualFluido_CalcTemp.Length; n++)
        {
            TactualPerdidas[n] = TactualFluido_CalcTemp[n];
        }

        Hl_local = CalculoPerdidasMetal(TactualPerdidas, Tambiente_CalcTemp);

        aux2 = CalculoCaracteristicasDelFluido(TactualFluido_CalcTemp, Caudal_CalcTemp);
        pf_local = new double[aux2.GetLength(1)];
        Cf_local = new double[aux2.GetLength(1)];
        Ht_local = new double[aux2.GetLength(1)];
        for(n = 0; n < aux2.GetLength(1); n++) // Quiero como limite del bucle el numero de columnas
        {   // Ojo que aqui los datos se organizan por filas, y no por columnas (como es el caso del fichero Datos.csv)
            pf_local[n] = aux2[0, n];
            Cf_local[n] = aux2[1, n];
            Ht_local[n] = aux2[2, n];
        }

        // Ajuste de caudal // Explicar por que
        if(Caudal_CalcTemp < 0.5)
        {
            aux2 = CalculoCaracteristicasDelFluido(TactualFluido_CalcTemp, Caudal_CalcTemp*10);
            Caudal_CalcTemp = 0.35;
            pf_local = new double[aux2.GetLength(1)];
            Cf_local = new double[aux2.GetLength(1)];
            Ht_local = new double[aux2.GetLength(1)];
            for(n = 0; n < aux2.GetLength(1); n++) // Quiero como limite del bucle el numero de columnas
            {   // Ojo que aqui los datos se organizan por filas, y no por columnas (como es el caso del fichero Datos.csv)
                pf_local[n] = aux2[0, n];
                Cf_local[n] = aux2[1, n];
                Ht_local[n] = aux2[2, n];
            }
        }
        // se pasa el caudal de m3/h a m3/s  dividiendo por 3600
        q = Caudal_CalcTemp / 3600;

        // Cálculo de la radiación efectiva, que llega a cada trozo del campo. Solo
        // llega radiación a las partes activas. Las partes pasivas tienen un
        // coeficiente de pérdidas mucho menor y la radiación incidente es nula

        Tfluido_CalcTemp        = new double[TactualFluido_CalcTemp.Length];
        Tfluido_CalcTemp[0]     = TactualFluido_CalcTemp[0];
        TemperaturasCalcTemp    = new double[2, TactualFluido_CalcTemp.Length];
        Tmetal_CalcTemp         = new double[TactualMetal_CalcTemp.Length];
        TfluidoSinCorregir      = new double[TactualFluido_CalcTemp.Length];

        vectorradiacion = new double[135];
        for(n = 0; n < 64; n++)
        {
            vectorradiacion[n] = I_CalcTemp*G*Eficienciamedia_CalcTemp*Factorsombra_CalcTemp*Porcentajetubo;
        }
        for(n = 64; n < 135; n++)
        {
            vectorradiacion[n] = 0;
        }

        // Paso1: Cálculo de la temperatura del metal
        for(n = 0; n < 64; n++)
        {
            Tmetal_CalcTemp[n] = TactualMetal_CalcTemp[n] + (tint_CalcTemp/(pm*Cm*Am))
                                *(vectorradiacion[n] - Hl_local[n]*G*(TactualMetal_CalcTemp[n] - Tambiente_CalcTemp)
                                - L*Ht_local[n]*(TactualMetal_CalcTemp[n] - TactualFluido_CalcTemp[n]));
            TemperaturasCalcTemp[0, n] = Tmetal_CalcTemp[n];
        }

        // Paso2: Cálculo de la temperatura del fluido en estado estacionario (q=0)
        for(n = 0; n < 64; n++)
        {
            TfluidoSinCorregir[n] = TactualFluido_CalcTemp[n]
                                    + ((L*Ht_local[n]*tint_CalcTemp)/(Af*pf_local[n]*Cf_local[n]))
                                    *(TactualMetal_CalcTemp[n] - TactualFluido_CalcTemp[n]);
        }
        TemperaturasCalcTemp[1, 0] = TfluidoSinCorregir[0];
        
        // Paso3: Corrección de la temperatura del fluido con la energía por transporte de caudal
        for(n = 1; n < 64; n++) // Explicar por que aqui el indice de inicio es 1 y no 0;
        {
            Tfluido_CalcTemp[n] = TfluidoSinCorregir[n] - (q*tint_CalcTemp)/(Af*incx[n-1])
                                    *(TfluidoSinCorregir[n] - TfluidoSinCorregir[n-1]);
            Debug.Log("Tfluido captador: "+Tfluido_CalcTemp[n]);
            TemperaturasCalcTemp[1, n] = Tfluido_CalcTemp[n];
        }
        

        // Tubería

        for(n = 64; n < 254; n++)
        {
            Tmetal_CalcTemp[n] = TactualMetal_CalcTemp[n];
            TemperaturasCalcTemp[0, n] = Tmetal_CalcTemp[n];
        }

        // Paso2: Cálculo de la temperatura del fluido en estado estacionario (q=0)
        for(n = 64; n < 254; n++)
        {
            TfluidoSinCorregir[n] = TactualFluido_CalcTemp[n] - L2*Hl_local[n]*(TactualFluido_CalcTemp[n] - Tambiente_CalcTemp)
                                    *tint_CalcTemp/(Af2*pf_local[n]*Cf_local[n]);
        }

        // Paso3: Corrección de la temperatura del fluido con la energía por transporte de caudal
        for(n = 64; n < 254; n++)
        {
            Tfluido_CalcTemp[n] = TfluidoSinCorregir[n] - (q*tint_CalcTemp)/(Af2*incx[n-1])
                                    *(TfluidoSinCorregir[n] - TfluidoSinCorregir[n-1]);
            TemperaturasCalcTemp[1, n] = Tfluido_CalcTemp[n];
        }

        return TemperaturasCalcTemp;
    }

    static void ObtenerDatos()
    {
        string[] FactorsombraCSV = CSVReader.LeerCSV(@"C:\Users\Pablo\OneDrive - UNIVERSIDAD DE SEVILLA\TFG\FSCmodel\FSCmodelUnity\Assets\DataFiles\Factorsombra.csv");
        string[] PorcentajetuboCSV = CSVReader.LeerCSV(@"C:\Users\Pablo\OneDrive - UNIVERSIDAD DE SEVILLA\TFG\FSCmodel\FSCmodelUnity\Assets\DataFiles\Porcentajetubo.csv");
        string[] TactualFluidoCSV = CSVReader.LeerCSV(@"C:\Users\Pablo\OneDrive - UNIVERSIDAD DE SEVILLA\TFG\FSCmodel\FSCmodelUnity\Assets\DataFiles\TactualFluido.csv");
        string[] TactualMetalCSV = CSVReader.LeerCSV(@"C:\Users\Pablo\OneDrive - UNIVERSIDAD DE SEVILLA\TFG\FSCmodel\FSCmodelUnity\Assets\DataFiles\TactualMetal.csv");

        bool success;
        int longitudBucle = FactorsombraCSV.GetLength(0);   // Importante que todos los ficheros de datos tengan la misma longitud
                                                            // o, al menos, que tome la longitud del más pequeño

        Factorsombra1       = new double[FactorsombraCSV.GetLength(0)];
        Porcentajetubo1     = new double[PorcentajetuboCSV.GetLength(0)];
        TactualFluido1      = new double[TactualFluidoCSV.GetLength(0)];
        TactualMetal1       = new double[TactualMetalCSV.GetLength(0)];

        for(int contador = 0; contador < longitudBucle; contador++)
        {
            success = Double.TryParse(FactorsombraCSV[contador], out Factorsombra1[contador]);      // h
            if(!success) Debug.Log("Error al convertir FactorsombraCSV a formato double, en el indice: " + contador);
            success = Double.TryParse(PorcentajetuboCSV[contador], out Porcentajetubo1[contador]);    // m3/s
            if(!success) Debug.Log("Error al convertir PorcentajetuboCSV a formato double, en el indice: " + contador);
            success = Double.TryParse(TactualFluidoCSV[contador], out TactualFluido1[contador]);  // ºC
            if(!success) Debug.Log("Error al convertir TactualFluidoCSV a formato double, en el indice: " + contador);
            success = Double.TryParse(TactualMetalCSV[contador], out TactualMetal1[contador]);   // ºC
            if(!success) Debug.Log("Error al convertir TactualMetalCSV a formato double, en el indice: " + contador);
        }
    }

}

