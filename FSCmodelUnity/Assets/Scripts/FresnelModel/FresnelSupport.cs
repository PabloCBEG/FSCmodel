using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class FresnelSupport
{
    // Aqui incluimos todas las funciones que tienen que ver con el modelo físico del captador Fresnel.

    public static double CalculoDiaJuliano(double anyo, double mes, double dia)
    {
        // Dias del mes para cada mes de un anyo normal
        int[] tabla = new int[] {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        // Dias del mes para cada mes de un anyo bisiesto
        int[] tablabisiesto = new int[] {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        int suma, i;
        double diajuliano;

        // Calculo del dia juliano segun la fecha proporcionada
        // Utilizamos calendario juliano por su continuidad en el tiempo frente al gregoriano

        // Comprobacion de si el anyo es bisiesto
        // En C# los arrays empiezan con indice 0
        if( anyo % 4 == 0 && ( anyo % 100 > 0 || anyo % 400 == 0 ))
        {
            /*See:
            https://learn.microsoft.com/en-us/office/troubleshoot/excel/determine-a-leap-year
            for an explanation of this method. */
            tabla[1] = tablabisiesto[1];
        }
        suma = 0;

        // Calculo del dia juliano
        for(i = 0; i < mes - 1; i++) // Ojo a ver si el mes es el correcto
        {
            suma = suma + tabla[i];
        }

        return diajuliano = suma + dia;
    }

    public static double CalculoAnguloDiario(double Diajuliano)
    {
        double anguloDiario;

        return anguloDiario = (2*Math.PI/365)*(Diajuliano - 1);
    }

    public static double CalculoHoraSolar(  double hora_HoraSolar,
                                            double minuto_HoraSolar,
                                            double mes_HoraSolar,
                                            double angulodiario_HoraSolar)
    {
        double Longitudmeridianocentralhuso, diferenciaporposicion, diferenciapormes, diferenciatotal;
        double horaactualoficial, Et, Ts;

        horaactualoficial = hora_HoraSolar;

        // Hora oficial en decimal
        horaactualoficial = horaactualoficial + minuto_HoraSolar/60;
        Longitudmeridianocentralhuso = 0;   // Meridiano de Greenwich

        // diferencia por posición en minutos
        diferenciaporposicion = Math.Abs((Fresnel.Longitud - Longitudmeridianocentralhuso)*4);

        // diferencia por horario de verano        
        if(mes_HoraSolar < 4 || mes_HoraSolar > 10)
            diferenciapormes = 60;  // 1 hora si no estamos entre abril y octubre
        else
            diferenciapormes = 120; // 2 horas entre abril y octubre

        // Ecuación del tiempo
        Et = (0.000075+0.001868*Math.Cos(angulodiario_HoraSolar)
            - 0.032077*Math.Sin(angulodiario_HoraSolar)
            - 0.014615*Math.Cos(2*angulodiario_HoraSolar)
            - 0.04089*Math.Sin(2*angulodiario_HoraSolar))*229.18; //En minutos

        // diferencia total en minutos
        diferenciatotal = diferenciapormes + diferenciaporposicion - Et;

        return Ts = horaactualoficial - (diferenciatotal/60);
    }

    public static double eficienciaGeoYSombras( double AnguloDiario_EficienciaGeoYSombras,
                                                double Ts_EficienciaGeoYSombras)
    {
        double[][]  resultadoVectorPosicionEspejos,
                    resultadoCalculoInclinacionEspejos,
                    resultadoCalculoEficienciaGeometrica;
        double[,]   r2d;
        double[]    alpha, theta, beta, kappa, inclinacion, factorsombra,
                    resultadoAzimut, eficienciageometrica, i2d, i3d, r2d0, r2d1;
        double      modi2d, gamma, Pps, alturasolar, azimut, Factorsombra;
        int         numeroFilasEspejos, j;

        // Definición de posición por fila espejos y otros datos físicos de la planta
        double[] XE = {-3.5, -2.8, -2.1, -1.4, -0.7, 0, 0.7, 1.4, 2.1, 2.8, 3.5};
        //private double Xtubo = 0;
        double Anchuraespejos = 0.25;
        int Ytubo = 4, Areatotal = 352, longitudPlanta = 64;

        numeroFilasEspejos = XE.Length;
        r2d0 = new double[numeroFilasEspejos];
        r2d1 = new double[numeroFilasEspejos];
        r2d = new double[2, numeroFilasEspejos];

        resultadoAzimut = CalculoAzimut(Ts_EficienciaGeoYSombras);
        alturasolar = resultadoAzimut[0];
        azimut      = resultadoAzimut[1];

        resultadoVectorPosicionEspejos = CalculoVectorPosicion2DEspejos(alturasolar, azimut, XE, Ytubo);
        i2d     = resultadoVectorPosicionEspejos[0];
        i3d     = resultadoVectorPosicionEspejos[1];
        r2d0    = resultadoVectorPosicionEspejos[2];
        r2d1    = resultadoVectorPosicionEspejos[3];
        modi2d  = resultadoVectorPosicionEspejos[4][0];
        gamma   = resultadoVectorPosicionEspejos[5][0];
        for(j = 0; j < r2d0.Length; j++)
        {
            r2d[0, j] = r2d0[j];
            r2d[1, j] = r2d1[j];
        }

        resultadoCalculoInclinacionEspejos = CalculoInclinacionEspejos(XE, i2d, r2d, modi2d, gamma);
        alpha       = resultadoCalculoInclinacionEspejos[0];
        theta       = resultadoCalculoInclinacionEspejos[1];
        beta        = resultadoCalculoInclinacionEspejos[2];
        kappa       = resultadoCalculoInclinacionEspejos[3];
        inclinacion = resultadoCalculoInclinacionEspejos[4];
        factorsombra = resultadoCalculoInclinacionEspejos[5];

        factorsombra = CalculoFactorSombra(XE, inclinacion, Anchuraespejos, gamma, factorsombra);

        resultadoCalculoEficienciaGeometrica = CalculoEficienciaGeometrica(factorsombra, Areatotal, inclinacion, kappa, gamma);
        Pps                     = resultadoCalculoEficienciaGeometrica[0][0];
        eficienciageometrica    = resultadoCalculoEficienciaGeometrica[1];

        Factorsombra = CalculoParteSombreada(inclinacion, alturasolar, i3d, longitudPlanta, Ytubo, eficienciageometrica, Pps, azimut, XE);

        return Factorsombra;
    }

    private static double[] CalculoAzimut(double Ts)
    {
        // Compmutes the azimuth and "solar height" for given data "solar time"
        double[]    resultado;
        double      angulodiario = Fresnel.angulodiario1;  // This is to be updated on every iteration (in fact, only when angulodiario changes)
        double      alturasolar, azimut, sinazimut, declinacion, anghorario, latitud;

        latitud = Fresnel.Latitud;

        // declinación calculada por Spencer 
        declinacion =   0.006918 - 0.399912*Math.Cos(angulodiario)
                        + 0.070257*Math.Sin(angulodiario)
                        - 0.006758*Math.Cos(2*angulodiario)
                        + 0.000907*Math.Sin(2*angulodiario)
                        - 0.002697*Math.Cos(3*angulodiario)
                        + 0.00148*Math.Sin(3*angulodiario);

        // ángulo horario                                        
        anghorario = ((Ts - 12)*15)*(Math.PI/180); // en radianes

        // cálculo de la altura y azimut solar                   
        alturasolar = Math.Asin(Math.Sin(declinacion)*Math.Sin((latitud*Math.PI)/180)
                    + Math.Cos(declinacion)*Math.Cos((latitud*Math.PI)/180)*Math.Cos(anghorario));
                
        azimut      = Math.Acos((Math.Sin(alturasolar)*Math.Sin((latitud*Math.PI)/180)
                    - Math.Sin(declinacion))/(Math.Cos(alturasolar)*Math.Cos((latitud*Math.PI)/180)));

        // Cálculo del signo del azimut solar                     
        sinazimut = (Math.Cos(declinacion)*Math.Sin(anghorario))/(Math.Cos(alturasolar));

        // si la cantidad es positiva el ángulo es negativo                
        if(sinazimut > 0) azimut = - azimut;

        return resultado = new double[] {alturasolar, azimut};
    }

    private static double[][] CalculoVectorPosicion2DEspejos(double alturasolar, double azimut, double[] XE, int Ytubo)
    {
        double[][]  resultado;
        double[]    i2d, i3d, r2d0, r2d1;
        double      oriplanta = Fresnel.oriplanta;
        double      gamma, modi2d;
        int         numeroFilasEspejos, j;

        numeroFilasEspejos  = XE.Length; // This is not changing (for this Solar plant)
        resultado           = new double[6][];
        r2d0                = new double[numeroFilasEspejos];
        r2d1                = new double[numeroFilasEspejos];
        i2d                 = new double[2];
        i3d                 = new double[3];

        // Modelo 2D, cálculo del vector i2d solar
        // No usamos vectores, sino arrays
        if(azimut > 0)
        {
            i2d = new double[] {Math.Cos(alturasolar)*Math.Cos(azimut + oriplanta), Math.Sin(alturasolar)};
            i3d = new double[] {i2d[0], i2d[1], Math.Cos(alturasolar)*Math.Sin(azimut + oriplanta)};
        }
        else if(azimut<0 && Math.Abs(azimut)<=oriplanta)
        {
            i2d = new double[] {Math.Cos(alturasolar)*Math.Cos(oriplanta - Math.Abs(azimut)), Math.Sin(alturasolar)};
            i3d = new double[] {i2d[0], i2d[1], Math.Cos(alturasolar)*Math.Sin(oriplanta - Math.Abs(azimut))};
        }
        else if(azimut<0 && Math.Abs(azimut)>oriplanta)
        {
            i2d = new double[] {Math.Cos(alturasolar)*Math.Cos(Math.Abs(azimut) - oriplanta), Math.Sin(alturasolar)};
            i3d = new double[] {i2d[0], i2d[1], Math.Cos(alturasolar)*Math.Sin(Math.Abs(azimut) - oriplanta)};
        }

        modi2d = Math.Sqrt(Math.Pow(i2d[0], 2) + Math.Pow(i2d[1], 2));
        
        for(j = 0; j < i3d.Length; j++)
        {
            i3d[j] = i3d[j]*(1/Math.Sqrt(Math.Pow(i3d[0], 2) + Math.Pow(i3d[1], 2) + Math.Pow(i3d[2], 2)));
        }

        gamma = Math.Atan(i2d[1]/i2d[0]);

        if(alturasolar <= 0) gamma = 0; //gamma=alturasolar;

        // Cálculo del vector r2d de posiciones de espejos.
        for(j = 0; j < numeroFilasEspejos; j++)
        {
            r2d0[j] = - XE[j]/Math.Sqrt(Math.Pow(XE[j], 2) + Math.Pow(Ytubo, 2));
            r2d1[j] = Ytubo/Math.Sqrt(Math.Pow(XE[j], 2) + Math.Pow(Ytubo, 2));
        }

        resultado[0] = i2d;
        resultado[1] = i3d;
        resultado[2] = r2d0; // r2d.CopyTo(resultado, 2);
        resultado[3] = r2d1;
        resultado[4] = new double[] {modi2d};
        resultado[5] = new double[] {gamma};

        return resultado;
    }

    private static double[][] CalculoInclinacionEspejos(double[] XE, double[] i2d, double[,] r2d, double modi2d, double gamma)
    {
        // Cálculo del ángulo que forma cada espejo con la horizontal beta.

        double[][]  resultado;
        double[]    alpha, theta, beta, kappa, inclinacion, factorsombra;
        int         numeroFilasEspejos = XE.Length, j;

        resultado   = new double[6][];
        alpha       = new double[numeroFilasEspejos]; theta         = new double[numeroFilasEspejos];
        beta        = new double[numeroFilasEspejos]; kappa         = new double[numeroFilasEspejos];
        inclinacion = new double[numeroFilasEspejos]; factorsombra  = new double[numeroFilasEspejos];

        for(j = 0; j < numeroFilasEspejos; j++)
        {
            beta[j] = Math.Atan(r2d[1, j]/Math.Abs(r2d[0, j]));
        }

        // Cálculo del ángulo alpha y theta
        for(j = 0; j < numeroFilasEspejos; j++)
        {
            alpha[j] = Math.Acos((r2d[0, j]*i2d[0] + r2d[1, j]*i2d[1])/modi2d);
            theta[j] = alpha[j]/2;
        }

        // Distintos casos, para cada espejo, primero para las 
        // primeras 5 filas. Cálculo de inclinación
        for(j = 0; j < beta.Length; j++)
        {
            inclinacion[j]  = 1;
            factorsombra[j] = 0;    // Inicializacion de factorsombra. Lo usamos mas adelante.
            kappa[j]        = 1;    // Inicializacion de kappa. Lo usamos mas adelante.
        }

        for(j = 0; j<5; j++)
        {
            if(gamma >= beta[j])
            {
                inclinacion[j] = Math.PI/2 - beta[j] - theta[j];
            }
            if(gamma < beta[j])
            {
                inclinacion[j] = Math.PI/2 - beta[j] + theta[j];
            }
        }
        inclinacion[5] = theta[5];
        for(j = 6; j < numeroFilasEspejos; j++)
        {
            inclinacion[j] = beta[j] + theta[j] - Math.PI/2;
        }

        resultado[0] = alpha;
        resultado[1] = theta;
        resultado[2] = beta;
        resultado[3] = kappa;
        resultado[4] = inclinacion;
        resultado[5] = factorsombra;

        return resultado;
    }

    private static double[] CalculoFactorSombra(double[] XE,
                                                double[] inclinacion,
                                                double Anchuraespejos,
                                                double gamma,
                                                double[] factorsombra)
    {
        // Una vez calculados los ángulos necesarios, se calculan
        // los factores de sombra y la eficiencia geométrica

        double[]    D, E;
        double      Pii, pj, dx, dy, modulo, w, a, b;
        int         j, k;

        for(k = 0; k < 10; k++)   // ¿Por que se hacen 10 iteraciones?
        {
            j = k + 1;
            Pii = inclinacion[k];
            pj = inclinacion[j];
            
            // analizamos los casos:
            // 1º Si ambas inclinaciones son positivas
            // 2º Si la primera positiva y la segunda negativa
            // 3º Si las dos son negativas
            if(Pii > 0 && pj > 0)
            {
                D = new double[] {XE[j] - Anchuraespejos*Math.Cos(pj), Anchuraespejos*Math.Sin(pj)};
                E = new double[] {XE[k] + Anchuraespejos*Math.Cos(Pii), -Anchuraespejos*Math.Sin(Pii)};
                
                // ángulo w
                dx  = Math.Abs(D[0] - E[0]);
                dy  = Math.Abs(D[1] - E[1]);
                modulo = Math.Sqrt(Math.Pow(dx, 2) + Math.Pow(dy, 2));
                w   = Math.Atan(dy/dx);
                a   = w - gamma;
                b   = Pii + gamma;

                if(gamma < w)
                {
                    factorsombra[k] = modulo*Math.Sin(a)/Math.Sin(b);
                }
            }
            if(Pii > 0 && pj < 0)
            {
                D = new double[] {XE[j] + Anchuraespejos*Math.Cos(pj), Anchuraespejos*Math.Sin(pj)};
                E = new double[] {XE[k] + Anchuraespejos*Math.Cos(Pii), -Anchuraespejos*Math.Sin(Pii)};
                
                // ángulo w
                dx  = Math.Abs(D[0] - E[0]);
                dy  = Math.Abs(D[1] - E[1]);
                modulo = Math.Sqrt(Math.Pow(dx, 2) + Math.Pow(dy, 2));
                w   = Math.Atan(dy/dx);
                a   = w - gamma;
                b   = Pii + gamma;

                if(gamma < w)
                {
                    factorsombra[k] = modulo*Math.Sin(a)/Math.Sin(b);
                }
            }
            if(Pii < 0 && pj < 0)
            {
                factorsombra[k] = 0;
            }
        }

        return factorsombra;
    }

    private static double[][] CalculoEficienciaGeometrica(  double[] factorsombra,
                                                            int Areatotal,
                                                            double[] inclinacion,
                                                            double[] kappa,
                                                            double gamma)
    {
        double[][]  resultado;
        double[]    eficienciageometrica;
        double      Pps;
        int         j;

        eficienciageometrica    = new double[kappa.Length];
        resultado               = new double[2][];

        for (j = 0; j < factorsombra.Length; j++)
        {
            factorsombra[j] = Math.Abs(factorsombra[j])*32;
        }
        Pps = 1 - factorsombra.Sum()/Areatotal;
        
        // Se halla ahora la parte del tubo que está
        // en sombra. Angulo kappa
        for(j = 0; j < 11; j++)
        {
            if(inclinacion[j] > 0)
            {
                kappa[j] = Math.PI/2 - gamma - inclinacion[j];
            }
            if(inclinacion[j] < 0)
            {
                kappa[j] = Math.PI/2 + inclinacion[j] - gamma;
            }
        }
        for (j = 0; j < kappa.Length; j++)
        {
            eficienciageometrica[j] = Math.Abs(Math.Cos(kappa[j]));
        }

        resultado[0] = new double[] {Pps};
        resultado[1] = eficienciageometrica;

        return resultado;
    }

    private static double CalculoParteSombreada(double[] inclinacion,
                                                double alturasolar,
                                                double[] i3d,
                                                int longitudPlanta,
                                                double Ytubo,
                                                double[] eficienciageometrica,
                                                double Pps, double azimut, double[] XE)
    {
        double[,]   ene, C;
        double[]    aux, ff, partesombra;
        double      Factorsombra, tubosombra, fulit, Eficienciatubosombra;
        int         ZR, j;
        int         numeroFilasEspejos = XE.Length;

        ene         = new double[3, numeroFilasEspejos];
        C           = new double[3, numeroFilasEspejos];
        ff          = new double[numeroFilasEspejos];
        partesombra = new double[numeroFilasEspejos];

        // Modelo 3D, parte del tubo que está en sombra
        for(j = 0; j < numeroFilasEspejos; j++) // Justificar por que 11 iteraciones aqui
        {
            ene[0, j] = Math.Sin(inclinacion[j]);
            ene[1, j] = Math.Cos(inclinacion[j]);
            ene[2, j] = 0;
        }

        if(azimut > 0) ZR = longitudPlanta;
        else ZR = 0;

        for(j = 0; j < numeroFilasEspejos; j++)
        {
            aux = new double[] {ene[0, j], ene[1, j], ene[2, j]};
            C[0, j]    = (double)CrossProduct3D(i3d, aux).GetValue(0);
            C[1, j]    = (double)CrossProduct3D(i3d, aux).GetValue(1);
            C[2, j]    = (double)CrossProduct3D(i3d, aux).GetValue(2);
            ff[j]   = Math.Pow((C[0, j]/Math.Cos(inclinacion[j])), 2); // Explicar por que C[0]
            partesombra[j] = ZR - Math.Sqrt((Math.Pow(XE[j], 2) + Math.Pow(Ytubo, 2))*ff[j]/(1 - ff[j]));
        }
        fulit = partesombra.Average();
        
        // Comprobamos si hay sombra efectivamente o no.
        if(ZR == 0)
        {
            tubosombra = Math.Abs(fulit/longitudPlanta);
        }
        if(ZR == longitudPlanta)
        {
            tubosombra = (longitudPlanta - fulit)/longitudPlanta;
        }
        // Hay que poner una alternativa, pues no es posible que tubosombra no tome ningun valor.
        // Vamos a considerar que, si no se da ninguno de los dos casos previos, tubosombra = 0;
        else tubosombra = 0;

        Eficienciatubosombra = 1 - tubosombra;

        if(Eficienciatubosombra > 1)
        {
            Eficienciatubosombra = 1;
        }
        if(Eficienciatubosombra < 0)
        {
            Eficienciatubosombra = 0;
        }

        Factorsombra = Eficienciatubosombra*eficienciageometrica.Max()*Pps;
        if(Factorsombra < 0 || alturasolar < 0)
        {
            Factorsombra = 0;
        }
        if(Factorsombra > 1)
        {
            Factorsombra = 1;
        }

        return Factorsombra;
    }

    private static double[] CrossProduct3D(double[] a, double[] b)
    {   // Performs cross-product for 3D arrays (not necessarily vectors, that's the point)
        double[] result = new double[3];
        result[0] = a[1] * b[2] - a[2] * b[1];
        result[1] = a[2] * b[0] - a[0] * b[2];
        result[2] = a[0] * b[1] - a[1] * b[0];
        return result;
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

        // Debug.Log("Temperatura fluido 33: "+TactualFluido_CalcTemp[33]);

        Af  = Fresnel.Af;  Af2  = Fresnel.Af2;  Am  = Fresnel.Am;  Am2 = Fresnel.Am2;
        G   = Fresnel.G;   L    = Fresnel.L;    L2  = Fresnel.L2;
        pm  = Fresnel.pm;  Cm   = Fresnel.Cm;
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
            // Debug.Log("Tfluido captador: "+Tfluido_CalcTemp[n]);
            TemperaturasCalcTemp[1, n] = Tfluido_CalcTemp[n];
        }
        

        // Tubería

        for(n = 64; n < 254; n++)
        {
            // En realidad la temperatura del metal en el sector de tubería fuera del captador
            // no nos importa, por eso asignamos aquí temperatura constante.
            // Tmetal_CalcTemp[n] = TactualMetal_CalcTemp[n];
            Tmetal_CalcTemp[n] = TactualMetal_CalcTemp[n] + (tint_CalcTemp/(pm*Cm*Am2))
                                *(vectorradiacion[n] - Hl_local[n]*Math.PI*0.16*(TactualMetal_CalcTemp[n] - Tambiente_CalcTemp)*0.1
                                - L2*0.01*Ht_local[n]*(TactualMetal_CalcTemp[n] - TactualFluido_CalcTemp[n]));
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
        double[] Hl;

        //Tambiente = GetComponent<CalcTemp>().Tambiente_CalcTemp;
        //tamanyoT_PerdidasMetal = GetComponent<CalcTemp>().TactualPerdidas.Length;
        //T = GetComponent<CalcTemp>().TactualPerdidas;
        tamanyoT_PerdidasMetal = T.Length;
        Hl = new double[tamanyoT_PerdidasMetal];

        Sup = 64*11*0.5;

        // Calculamos por metro cuadrado de espejo para perdidas de tubo             
        for(m = 0; m < tamanyoT_PerdidasMetal; m++)
        {
            Hl[m] = (0.45659247191149893/Sup)*(T[m] - Tambiente) - 1.062045206593238/Sup;
        }

        // Metros 0-64 de la tuberia se corresponden al captador;
        // Metros 65 en adelante, al resto de la tuberia, intercambiador de calor...
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

    
}