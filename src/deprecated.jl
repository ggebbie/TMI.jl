
# """
#     function adjustboundarycondition(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where N1, N2, T <: Real

#     adjust all boundary conditions b that are described in u
# """
# function adjustboundarycondition(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where {N1, N2, T <: Real}

#     bnew = deepcopy(b)
#     ukeys = keys(u)
#     for ukey in keys(u)
#         bnew[ukey].tracer[bnew[ukey].wet] += u[ukey].tracer[bnew[ukey].wet] 
#     end
#     return bnew
# end

# function unvec!(u::BoundaryCondition{T},uvec::Vector{T}) where T <: Real
#     I = findall(u.wet)
#     counter = 0
#     for i in I
#         counter +=1
#     # doesn't work to pass u back
#         #u.tracer[u.wet] = uvec
#         u.tracer[i] = uvec[counter]
#     end
# end
# function unvec!(u::NamedTuple,uvec::Vector) #where {N, T <: Real}

#     #counter = 0
#     for v in u
#         # if v isa BoundaryCondition
#         #     wet = v.wet
#         # elseif v isa Field
#         #     wet = v.γ.wet
#         # end
#         n = sum(wet(v))
#         v.tracer[wet(v)] = uvec[counter+1:counter+n]
#         counter += n
#     end
# end


# function zerotemplate(utemplate,uvec)
#     tracer = zeros(wet(utemplate))
#     tracer[wet(utemplate)] = uvec
#     return tracer
# end

# """
#     function unvec(u,uvec)

#     Replace u with new u
#     Undo the operations by vec(u)
#     Needs to update u because attributes of 
#     u need to be known at runtime.
# """
# function unvec(utemplate::Union{Field{T},BoundaryCondition{T}},uvec::Vector{T}) where T <: Real
#     tracer = zerotemplate(utemplate,uvec)
#     u = BoundaryCondition(tracer,utemplate.i,utemplate.j,utemplate.k,utemplate.dim,utemplate.dimval,wet(utemplate))
#     return u
# end
# function unvec(utemplate::BoundaryCondition{T},uvec::Vector{T}) where T <: Real
#     tracer = zeros(utemplate.wet)
#     tracer[utemplate.wet] = uvec
#     u = BoundaryCondition(tracer,utemplate.i,utemplate.j,utemplate.k,utemplate.dim,utemplate.dimval,utemplate.wet)
#     return u
# end
# function unvec(utemplate::Field{T},uvec::Vector{T}) where T <: Real
#     tracer = zeros(utemplate.γ)
#     tracer.tracer[utemplate.γ.wet] .= uvec
#     #u = BoundaryCondition(tracer,utemplate.i,utemplate.j,utemplate.k,utemplate.dim,utemplate.dimval,utemplate.wet)
#     return tracer
# end
# function unvec(utemplate::NamedTuple{<:Any,NTuple{N,BoundaryCondition{T}}},uvec::Vector{T}) where {N, T <: Real}
#     counter = 0
#     vals = Vector{BoundaryCondition}(undef,N)
#     for (ii,vv) in enumerate(utemplate)
#         n = sum(vv.wet)
#         vals[ii] = unvec(vv, uvec[counter+1:counter+n])
#         counter += n
#     end
#     u = (;zip(keys(utemplate), vals)...)
#     return u
# end
# function unvec(utemplate::NamedTuple{<:Any,NTuple{N,Field{T}}},uvec::Vector{T}) where {N, T <: Real}
#     counter = 0
#     vals = Vector{Field}(undef,N)
#     for (ii,vv) in enumerate(utemplate)
#         n = sum(vv.γ.wet)
#         vals[ii] = unvec(vv, uvec[counter+1:counter+n])
#         counter += n
#     end
#     u = (;zip(keys(utemplate), vals)...)
#     return u
# end
# function unvec(utemplate,uvec::Vector{T}) where {T <: Real}
#     counter = 0
#     vals = Vector{Any}(undef,length(utemplate))
#     for (ii,vv) in enumerate(utemplate)
#         if vv isa Field
#             n = sum(vv.γ.wet)
#         elseif vv isa BoundaryCondition
#             n = sum(vv.wet)
#         end
#         vals[ii] = unvec(vv, uvec[counter+1:counter+n])
#         counter += n
#     end
#     u = (;zip(keys(utemplate), vals)...)
#     return u
# end

