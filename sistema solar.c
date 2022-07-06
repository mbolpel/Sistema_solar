/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// SIMULACIÓN DEL SISTEMA SOLAR                                                //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

//BIBLIOTECAS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//CONSTANTES DEL PROGRAMA
#define h 0.001 //Intervalo de tiempo
#define N 2000000 //Numero cálculos del algoritmo

int main()
{
    //DEFINICIÓN DE VARIABLES
    double m, x, y, vx, vy; //masa, posición en eje x e y, velocidad en eje x e y
    double Et, Ec, Ep;      //Energías total, cinética y potencial
    int n;                  //contador
    double M[9];
    double pos[9][2], vel[9][2], acel[9][2];     //Matrices de datos (psoición, velocidad y aceleración)
    double w[9][2];        //Matriz del algoritmo de Verlet 
    double aux[9][2];      //Matriz auxiliar que guarda el valor de la posición de la iteración anterior
    double E[9];           //Vector de energía 
    FILE *f1,*f2,*f3;      //Ficheros

    //PASO 1: TOMA DE DATOS
    //1.1 Lectura de datos de los ficheros
    f1=fopen("distancias_ss.dat","r"); // en km
    f2=fopen("velocidades_ss.dat","r"); // en km/s
    f3=fopen("masas_ss.dat","r"); //en kg

    //1.2 Distancias (pasar los datos del fichero a la matriz de posiciones, cargando las distancias al sol en el eje de la x)
    n=0;
    while(fscanf(f1,"%lf",&x)!=EOF)
    {
        pos[n][0]=x;
        n++;
    }
    n=0;
    for(y=0; y<9; y++)
    {
        pos[n][1]=0;
        n++;
    }

    //1.3 Velocidad orbital (La velocidad es perpendicular a la posición, luego se cargan los datos en la componente y del vector)
    n=0;
    for(vx=0; vx<9; vx++)
    {
        vel[n][0]=0;
        n++;
    }
    n=0;
    while(fscanf(f2,"%lf",&vy)!=EOF)
    {
        vel[n][1]=vy;
        n++;
    }

    //1.4 Masas 
    n=0;
    while(fscanf(f3,"%lf",&m)!=EOF)
    {
        M[n]=m;
        n++;
    }

    //1.4 Mostrar los valores iniciales por pantalla----> Me sirve para comprobar que se han pasado bien los datos del fichero
    //Valores de la matriz de distancias(x,y)
    printf("\n\nMATRIZ DE DISTANCIAS\n\n");     

    for(int i=0;i<9;i++)                    
    {
        for(int j=0;j<2;j++)
        {
            printf("%lf\t",pos[i][j]);
        }
        printf("\n");
    }
    
    //Valores de la matriz de velocidades(x,y)
    printf("\n\nMATRIZ DE VELOCIDADES\n\n");        

    for(int i=0;i<9;i++)                    
    {
        for(int j=0;j<2;j++)
        {
            printf("%lf\t",vel[i][j]);
        }
        printf("\n");
    }

    //Valores del vector de masas
    /*printf("VECTOR DE MASAS\n\n");

    for(int i=0;i<9;i++)
    {
        printf("v=%lf\n";M[i]);
    }*/

    //Cerrar los ficheros
    fclose(f1);
    fclose(f2);
    fclose(f3);


    //PASO 2: ALGORITMO DE VERLET/////////////////////////////////////////////////////////////////////////////////////////////////
    // En este algoritmo primero necesito los valores iniciales de las posiciones, velocidades, aceleraciones y masas,          //
    // por ello necesito el cálculo del vector de aceleraciones iniciales (paso 2.1) y luego haré el cálculo de los distintos   //
    // valores t+h donde h es la variable  que indica el intervalo de tiempo, definida al principio del programa.               //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //2.1 Cálculo de la aceleración inicial
    // La i señala el planeta del que se está calculando la aceleración
    // La j señala el planeta respecto al que estamos viendo las distancias
    // La k señala a la coordenada (x o y) de las posiciones y de la aceleración
    for(int i=0;i<9;i++)
    {
        for(int j=0;j<9;j++)
        {
            if(j!=i)        //No se considera la interacción de un planeta consigo mismo
            {
                for(int k=0;k<2;k++)
                {
                    double distanciaplanetas;   //Para calcular la distancia entre dos planetas. 
                    double dist[2];             //Vector que guarda las coordenadas de la distancia para el eje x y el eje y
                                                
                    for(int a=0;a<2;a++) dist[a]=pos[i][a]-pos[j][a];

                    distanciaplanetas=sqrt(pow(dist[0],2)+pow(dist[1],2)); //La distancia entre planetas se calcula con la raiz de los cuadrados de sus coordenadas
                    //Mostrar las distancias entre planetas por pantalla
                    printf("\n\nDISTANCIAS DE LOS PLANETAS:\n\n"); 
                    printf("%lf\t",distanciaplanetas);
                    
                    acel[i][k]=acel[i][k]-M[j]*(pos[i][k]-pos[j][k])/(pow(distanciaplanetas,3));

                }
            }
        }
    }

    //Mostrar los valores de la matriz de aceleraciones(x,y)
    printf("\n\nMATRIZ DE ACELERACIONES\n\n");        

    for(int i=0;i<9;i++)                    
    {
        for(int j=0;j<2;j++)
        {
            printf("%lf\t",acel[i][j]);
        }
        printf("\n");
    }

    //2.2 Definición de ficheros donde se guardarán los datos
    FILE *f4, *f5, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f14;

    //Abrir los ficheros 
    f4=fopen("mercurio.dat","w");
    f5=fopen("venus.dat","w");
    f6=fopen("tierra.dat","w");
    f7=fopen("marte.dat","w");
    f8=fopen("jupiter.dat","w");
    f9=fopen("saturno.dat","w");
    f10=fopen("urano.dat","w");
    f11=fopen("neptuno.dat","w");
    f12=fopen("pos.dat","w");
    f13=fopen("energia.dat","w");
    f14=fopen("energiat.dat","w");

    int p;//Necesito un contador nuevo (p coge los valores que se obtienen en una de cada 1000 interacciones del bucle y las escribe en el fichero)
    p=0;

    //Primer nivel del bucle: contador de interacciones (k)
    for(int k=0;k<N;k++)
    {
        /*printf("\nValor de k=");   //Muestro la interacción en la que me encuentro
        printf("%i",k);

        Et=0.0;   //Pongo la energía total a cero
        */

        //Segundo nivel: i indica los planetas
        for (int i=0;i<9;i++)
        {
            //Tercer nivel: j varía la coordenada
            for(int j=0;j<2;j++)
            {
                pos[i][j]=pos[i][j]+h*vel[i][j]+pow(h,2)/2*acel[i][j];   //Cálculo de r(t+h)
                
                //2.3 Cálculo y guardado del periodo en su archivo correspondiente
                if(i>0 && pos[i][0]>0 && (pos[i][1]/aux[i][1])<=0)
                {
                    /*printf("\nEstamos en el bucle cuando k=");
                    printf("%i",k);
                    */

                    if(i==1) fprintf(f4,"%lf dias\n",k*(1.395/24)); //Mercurio
                    if(i==2) fprintf(f5,"%lf dias\n",k*(1.395/24)); //Venus
                    if(i==3) fprintf(f6,"%lf dias\n",k*(1.395/24)); //Tierra
                    if(i==4) fprintf(f7,"%lf dias\n",k*(1.395/24)); //Marte
                    if(i==5) fprintf(f8,"%lf dias\n",k*(1.395/24)); //Jupiter
                    if(i==6) fprintf(f9,"%lf dias\n",k*(1.395/24)); //Saturno
                    if(i==7) fprintf(f10,"%lf dias\n",k*(1.395/24)); //Urano
                    if(i==8) fprintf(f11,"%lf dias\n",k*(1.395/24)); //Neptuno
                }
                aux[i][j]=pos[i][j];                          //Guardo el valor de la posición en esta interacción
                //2.4 Guardado de las posiciones en una de cada 1000 interacciones
                if(p==1000) fprintf(f12,"lf\t",pos[i][j]);
                w[i][j]=vel[i][j]+h/2*acel[i][j];            //Cálculo de la velocidad w de Verlet
            }
        }

        //2.5 Calcular los nuevos valores de la aceleración
        //Iniciaizar las aceleraciones a cero
        for(int i=0;i<9;i++)
        {
            for(int j=0;j<2;j++) acel[i][j]=0.0;
        }
        
        //Mismo formato que el apartado 2.1
        for(int i=0;i<9;i++)
        {
            for(int j=0;j<9;j++)
            {
                if(j!=i)
                {
                    for(int a=0;a<2;a++)
                    {
                        double distanciaplanetas;
                        double dist[2];

                        for(int b=0;b<2;b++) dist[b]=pos[i][b]-pos[j][b];

                        distanciaplanetas=sqrt(pow(dist[0],2)+pow(dist[1],2));

                        acel[i][a]=acel[i][a]-M[j]*(pos[i][a]-pos[j][a])/pow(distanciaplanetas,3);
                    }
                }
            }
        }

        //2.6 Calcular las nuevas velocidades
        for(int i=0;i<9;i++)
        {
            for(int j=0;j<2;j++)
            {
                vel[i][j]=w[i][j]+h/2*acel[i][j];    
            }
        }

        //2.7 Cálculo de la energía
        //2.7.1 Guardo en el fichero de energías el valor de la interacción en la que me encuentro
        if(p==1000) fprintf(f13,"%i\t",k);

        //2.7.2 Cálculo de la energía total (sistema)
        for(int i=0;i<9;i++)
        {
            //Energía cinética
            double v;   //Defino lo que será el módulo de la velocidad
            v=sqrt(pow(vel[i][0],2)+pow(vel[i][1],2));

            Ec=1/2*M[i]*pow(v,2);

            //Energía potencial (como la energía potencial es una energía de interacción entre planetas se calcula la distancia entre planetas)
            Ep=0.0;

            double distanciaplanetas;
            double dist[2];

            for(int j=0;j<9;j++)
            {
                if(i!=j)
                {
                    for(int a=0;a<2;a++) dist[a]=pos[i][a]-pos[j][a];
                    distanciaplanetas=sqrt(pow(dist[0],2)+pow(dist[1],2));

                    Ep=Ep-M[i]*M[j]/distanciaplanetas;
                }
            }

            //Guardo la suma de energía potencial y cinética en el vector de energías
            E[i]=Ec+Ep;

            //Imprimo el valor de la energía de cada planeta cada 1000 iteracciones de k
            if(p==1000) fprintf(f13,"%E\t",E[i]);

            //Energía total del sistema
            Et=Et+E[i];
        }
        //Imprimo la energía total del sistema en el fichero
        if(p=1000)
        {
            fprintf(f14,"%i\t%E\n",k,Et);
            fprintf(f13,"\n");
            fprintf(f12,"\n");
            p=0;
        }
        p=p+1;
    }

    //Cierro todos los ficheros
    fclose(f4);
    fclose(f5);
    fclose(f6);
    fclose(f7);
    fclose(f8);
    fclose(f9);
    fclose(f10);
    fclose(f11);
    fclose(f12);
    fclose(f13);
    fclose(f14);

    return 0;
}