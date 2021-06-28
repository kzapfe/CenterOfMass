using JLD

## Auxiliary Functions
include("FuncionesTrayectorias.jl")


# Argument Filename of type jld is mandatory
numarg=length(ARGS)
if numarg<1
    error(" Argument Filename of type jld is mandatory")
else
    nombre=ARGS[1] 
end

nombregeneral=nombre[1:end-4]

nombresalida=string(nombregeneral, "-Tray.jld")

println("This shall be the outname: ", nombresalida)

if numarg==1
    inicadena=0
    fincadena=400
elseif numarg==2
    inicadena=0
    fincadena=ARGS[2]
elseif numarg==3
    inicadena=ARGS[2]
    fincadena=ARGS[3]
end

    


archivo=load(nombre)
DatosCMP=archivo["CMP"]
DatosCMN=archivo["CMN"]

tolerancia=3

CatenarioPositivo=encuentraTrayectorias(DatosCMP,tolerancia,inicadena,fincadena);
CatenarioNegativo=encuentraTrayectorias(DatosCMN,tolerancia,inicadena,fincadena);

println("Writing Output...")
save(nombresalida, "CatenarioNegativo", CatenarioNegativo, "CatenarioPositivo", CatenarioPositivo)
println("¡cha chan!")
