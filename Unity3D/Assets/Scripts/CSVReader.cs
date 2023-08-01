using UnityEngine;
using System;
using System.Linq;
using System.IO;
using System.Collections;
using System.Collections.Generic;

public static class CSVReader
{
    public static string[,] LeerCSV(string filePath)
    {
        int i = 0, j;       // Contadores para recorrer el array donde vamos a almacenar los valores.
        string[] values;    // Valores de las cadenas de cada línea del archivo .csv.
        string[,] valor;    // Array 2D de strings que contiene los valores del .csv organizados en matriz.
        
        int numberOfLines = File.ReadLines(filePath).Count();       //  Primero calculamos el número de líneas del archivo. En el caso
                                                                    // de nuestro proyecto, eso equivale al número de muestras reales.
        var firstLine = File.ReadLines(filePath).FirstOrDefault();  //  Ahora identificamos la primera línea.
        var numberOfColumns = firstLine.Split(',').Length;          // Con el objetivo de calcular el número de columnas que tiene el
                                                                    // archivo. Es decir, el número de variables medidas, en
                                                                    // nuestro caso.
        valor = new string[numberOfLines, numberOfColumns];         //  Inicializamos las dimensiones de nuestra variable "valor".

        if (File.Exists(filePath))
        {
            using (StreamReader reader = new StreamReader(filePath))
            {
                while (!reader.EndOfStream)             //  Mientras queden datos por leer:
                {
                    string line = reader.ReadLine();    //  Leemos la línea en curso,
                    values = line.Split(',');           // separamos los valores que contiene según las comas (importante que estén puestas),
                    j = 0;                              // reseteamos j a 0 antes de entrar en el bucle de asignación.
                    foreach (string value in values)    //  Como hemos separado los valores por las comas, ahora tenemos tantos strings
                    {                                   // como columnas tenga el archivo.
                        valor[i, j] = value;            //  Así que, cogemos cada string y lo asignamos a una posición del array 2D.
                        j++;
                    }
                    i++;
                }
            }
        }
        else
        {
            Debug.LogError("File not found: " + filePath);  //  Si no existe el archivo cuyo path nos han pasado, sacamos error por pantalla.
        }

        return valor;
    }
}

