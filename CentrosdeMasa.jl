module CentrosdeMasa


#=
Solo el calculo del centro de masa a partir
de un JLD que contenga un arreglo 3D que
tenga el CSD calculado.

Muy importante que los electrodos saturados
esten en cero.
=#


# Funciones

export vecindad8, ComponentesSP, ObtenComponentesyCM

function vecindad8(punto::Array)
    j=punto[1]
    k=punto[2]
    result=Set{Array{Int64,1}}()
    push!(result, [j-1,k-1])
    push!(result, [j-1,k])
    push!(result, [j-1,k+1])
    push!(result, [j,k-1])
    push!(result, [j,k+1])
    push!(result, [j+1,k-1])
    push!(result, [j+1,k])
    push!(result, [j+1,k+1])
    return result
end

function ComponentesSP(DatosSignados::Array)
    #Single pass method para sacar componentes disjuntos.
    lista=copy(DatosSignados)
    componentes=Set{Any}()
    while(length(lista)!=0)
        
        x=pop!(lista) #arranca el ULTIMO elemento de la lista
        listaprofundeza=Array{Int64}[]
        componentecurlab=Array{Int64}[]
        push!(listaprofundeza, x) #Pone elementos al FINAL de la lista
        push!(componentecurlab, x)    
        profundidad=0
        
        while ((length(listaprofundeza)!=0) && profundidad<1000)

            y=pop!(listaprofundeza)
            for v in vecindad8(y)
                if in(v, lista)
                    deleteat!(lista, indexin(Any[v], lista))
                    push!(listaprofundeza, v) 
                    profundidad+=1
                    push!(componentecurlab, v)
                end
            end
        end
        push!(componentes, componentecurlab)    
    end
        return componentes
    end
    



    
function ObtenComponentesyCM(Datos::Array, tini=1,tfini=tmax, epsilon=1.0)
    #CSD ahora no tiene orillas. Asi que toca adaptarse.
    (alto,ancho,lu)=size(Datos)
    #la cantidad minima de pixeles que tiene que tener un componente para
    #que lo tomemeos en cuenta
    tamano=3
    #Esto va a a ser el resultado de la funcion!
    #La llave es t
    #El contenido es la lista de CM.
    CMPositivo=Dict{Int, Array}()
    CMNegativo=Dict{Int, Array}()
    #Aqui empieza el circo
    for t=tini:tfini
        #iniciar variables vacias
        ActividadNegativa=Array{Int16}[]
        ActividadPositiva=Array{Int16}[]
        SpikeCountPositivo=zeros(alto,ancho)
        SpikeCountNegativo=zeros(alto,ancho)
        #Separamos pixeles positivos y negativos
        for j=1:alto,k=1:ancho
            if(Datos[j,k,t]<-epsilon)
                
                push!(ActividadNegativa, [j, k])
                SpikeCountNegativo[j,k]+=1
                
            elseif(Datos[j,k,t]>epsilon)
                
                push!(ActividadPositiva, [j, k])
                SpikeCountPositivo[j,k]+=1
                
            end
                
            end
            
            #Primero Negativo
        componentesneg=ComponentesSP(ActividadNegativa)
        centrosdemasaneg=[[0 0 0];]
            #=
        componentesneg/pos son
            conjuntos con las listas de elemenentos de los
        componentes en un instante dado. Se tiene que "cerar" siempre. 
       =#
        
        for p in componentesneg
            mu=length(p)        
            if mu>tamano
                masa=0.00
                x=0.00
                y=0.00
                for q in p
                    j=q[1]
                    k=q[2]
                    masalocal=Datos[j,k,t]
                    masa+=masalocal
                    x+=k*masalocal
                    y+=j*masalocal
                end
                x/=masa
                y/=masa
                A=[x y masa]
                centrosdemasaneg=vcat(centrosdemasaneg, A)
            end
        end
            centrosdemasaneg=centrosdemasaneg[2:end,:]
      #      println(t,"hola")
            CMNegativo[t]=centrosdemasaneg
 
        
        ##### Ahora lo posittivo (fuentes)
            componentespos=ComponentesSP(ActividadPositiva)               
            centrosdemasapos=[[0 0 0];]
            for p in componentespos
                mu=length(p)
                if mu>tamano
                    masa=0.00
                    x=0.00
                    y=0.00
                    for q in p
                        j=q[1]
                        k=q[2]
                        masalocal=Datos[j,k,t]
                        masa+=masalocal
                        x+=k*masalocal
                        y+=j*masalocal
                    end
                    x/=masa
                    y/=masa
                    A=[x y masa]
                    centrosdemasapos=vcat(centrosdemasapos, A)
                end
            end
        
        centrosdemasapos=centrosdemasapos[2:end,:]       
        CMPositivo[t]=centrosdemasapos
        
        end
    
    return (CMPositivo, CMNegativo)
    end
    
       

end # modulo
