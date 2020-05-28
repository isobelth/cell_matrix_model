parameter_dict = Dict("centre_1_position_x" => 1, "centre_1_position_y" => 2,
                      "centre_2_position_x" => 3, "centre_2_position_y" => 4,
                      "current_radius" => 5, "current_orientation" => 6, "current_deformation" => 7 )

maintained_parameters_dict = Dict("R_max" => 1, "dR_max" => 2, "dr_max" => 3,
                                  "dtheta_max" => 4, "dd_max" => 5, "initial_deformation_proportion" => 6)

cluster_parameters_dict = Dict("matrix_stiffness" => 1, "attractor_start_x" => 2, "attractor_start_y" => 3,
                                "attractor_end_x" => 4, "attractor_end_y" => 5, "attractor_strength" => 6)
function create_maintained_parameters(R_max::Float64, dR_max::Float64, dr_max::Float64, dtheta_max::Float64, dd_max::Float64, initial_deformation_proportion::Float64)
  maintained_parameters::Array{Float64, 1} = [R_max, dR_max, dr_max, dtheta_max, dd_max, initial_deformation_proportion]
  return maintained_parameters
end
function create_cluster_parameters(matrix_stiffness::Float64, attractor_start_x::Float64, attractor_start_y::Float64, attractor_end_x::Float64, attractor_end_y::Float64, attractor_strength::Float64)
    cluster_parameters::Array{Float64, 1} = [matrix_stiffness, attractor_start_x, attractor_start_y, attractor_end_x, attractor_end_y, attractor_strength]
    return cluster_parameters
end
#############FUNCTIONS FOR I CELLS###################
function create_i_cell(centre_1_position_x::Float64, centre_1_position_y::Float64, current_radius::Float64, current_orientation::Float64)
    i_cell::Array{Float64, 1} = [centre_1_position_x, centre_1_position_y, centre_1_position_x, centre_1_position_y, current_radius, current_orientation, 0.0]
    return i_cell
end
function migrate_i_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    direction_to_move = 2*pi*rand()
    amount_to_move = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    i_cell[parameter_dict["centre_1_position_x"]] += amount_to_move*cos(direction_to_move)
    i_cell[parameter_dict["centre_1_position_y"]] += amount_to_move*sin(direction_to_move)
    i_cell[parameter_dict["centre_2_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]]
    i_cell[parameter_dict["centre_2_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]]
    return i_cell
end
function grow_i_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    i_cell[parameter_dict["current_radius"]] += maintained_parameters[maintained_parameters_dict["dR_max"]]*rand()
    if i_cell[parameter_dict["current_radius"]] >  maintained_parameters[maintained_parameters_dict["R_max"]]
        i_cell[parameter_dict["current_radius"]] = maintained_parameters[maintained_parameters_dict["R_max"]]
    end
    return i_cell
end
###################FUNCTIONS FOR M CELLS#############
function create_m_cell(centre_1_position_x::Float64, centre_1_position_y::Float64, centre_2_position_x::Float64, centre_2_position_y::Float64, current_radius::Float64, current_orientation::Float64, current_deformation::Float64)
    m_cell::Array{Float64, 1} = [centre_1_position_x, centre_1_position_y, centre_2_position_x, centre_2_position_y, current_radius, current_orientation, current_deformation]
    return m_cell
end
function migrate_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    direction_to_move = 2*pi*rand()
    amount_to_move = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    m_cell[parameter_dict["centre_1_position_x"]] += amount_to_move*cos(direction_to_move)
    m_cell[parameter_dict["centre_1_position_y"]] += amount_to_move*sin(direction_to_move)
    m_cell[parameter_dict["centre_2_position_x"]] += amount_to_move*cos(direction_to_move)
    m_cell[parameter_dict["centre_2_position_y"]] += amount_to_move*sin(direction_to_move)
    return m_cell
end
function deform_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    amount_to_deform = maintained_parameters[maintained_parameters_dict["dd_max"]]*rand()
    m_cell[parameter_dict["current_deformation"]] += amount_to_deform
    d::Float64=m_cell[parameter_dict["current_deformation"]]
    R::Float64=m_cell[parameter_dict["current_radius"]]
    x::Float64=d/(2*R)
    max_area::Float64 = round(pi * maintained_parameters[maintained_parameters_dict["R_max"]] * maintained_parameters[maintained_parameters_dict["R_max"]], digits=4)
    if x<=1
        new_area::Float64 = round(2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d)), digits=3)
    else
        new_area::Inf16
    end
    while (abs(new_area-max_area))>0.5
        if new_area - max_area >0.5
            R-=0.001

            new_area = round(2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d)), digits=4)
            print("decrease radius, r=", R," new area=", new_area, "max area=", max_area, "\n")
        else
            R+=0.001
            new_area = round(2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d)), digits=4)
            print("increase radius, r=", R ," new area=", new_area,"max area=", max_area, "\n")
        end
    end
    m_cell[parameter_dict["current_radius"]] = round(R, digits=4)
    m_cell[parameter_dict["centre_1_position_x"]] += (amount_to_deform/2)*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_1_position_y"]] += (amount_to_deform/2)*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_x"]] += (amount_to_deform/2)*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_y"]] += (amount_to_deform/2)*sin(m_cell[parameter_dict["current_orientation"]])
    return m_cell
end
function rotate_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    amount_to_rotate::Float64 = maintained_parameters[maintained_parameters_dict["dtheta_max"]]*rand()
    m_cell[parameter_dict["current_orientation"]] += amount_to_rotate

    m_cell[parameter_dict["centre_1_position_x"]] += m_cell[parameter_dict["current_deformation"]]*0.5*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_1_position_y"]] += m_cell[parameter_dict["current_deformation"]]*0.5*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_x"]] -= m_cell[parameter_dict["current_deformation"]]*0.5*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_y"]] -= m_cell[parameter_dict["current_deformation"]]*0.5*sin(m_cell[parameter_dict["current_orientation"]])
    return m_cell
end
#############FUNCTIONS TO CHANGE CELL TYPE#######
function replace_i_cell_with_m_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    new_m_cell = create_m_cell(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    current_orientation = i_cell[parameter_dict["current_orientation"]]#####can update to divide in preferential directions
    R::Float64 = maintained_parameters[maintained_parameters_dict["R_max"]] ##max radius
    max_area = pi * R* R
    d::Float64 = maintained_parameters[maintained_parameters_dict["initial_deformation_proportion"]]*rand() ##initial displacement

    new_area::Float64 = round(2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d)), digits=3)
    while (abs(new_area-max_area))>0.5
        if new_area - max_area >0.1
            R-=0.01
            new_area = 2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d))
        else
            R+=0.01
            new_area = 2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d))
        end
    end
    new_m_cell[parameter_dict["current_radius"]] = round(R, digits=3)
    new_m_cell[parameter_dict["current_orientation"]] = current_orientation
    new_m_cell[parameter_dict["current_deformation"]] = d
    new_m_cell[parameter_dict["centre_1_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]] + (d)*0.5*cos(current_orientation)
    new_m_cell[parameter_dict["centre_1_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]] + (d)*0.5*sin(current_orientation)
    new_m_cell[parameter_dict["centre_2_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]] - (d)*0.5*cos(current_orientation)
    new_m_cell[parameter_dict["centre_2_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]] - (d)*0.5*sin(current_orientation)

    return new_m_cell
end
function replace_m_cell_with_two_i_cells(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    R = m_cell[parameter_dict["current_radius"]]
    current_orientation = m_cell[parameter_dict["current_orientation"]]
    new_i_cell_1_x_position = m_cell[parameter_dict["centre_1_position_x"]]
    new_i_cell_1_y_position = m_cell[parameter_dict["centre_1_position_y"]]
    new_i_cell_2_x_position = m_cell[parameter_dict["centre_2_position_x"]]
    new_i_cell_2_y_position = m_cell[parameter_dict["centre_2_position_y"]]
    new_i_cell_1::Array{Float64, 1} = create_i_cell(new_i_cell_1_x_position,new_i_cell_1_y_position,R,current_orientation)
    new_i_cell_2::Array{Float64, 1} = create_i_cell(new_i_cell_2_x_position,new_i_cell_2_y_position,R,current_orientation)
    new_i_cells::Array{Float64} = hcat(new_i_cell_1 , new_i_cell_2)
    #print("new i cells size", size(new_i_cells))
    return new_i_cells
end
############FUNCTION TO SET UP CLUSTER########
function initialise_cluster(i_cell::Array{Float64,1})
    cluster::Array{Float64, 1} = i_cell
    return cluster
end
function make_a_move(cluster::Array{Float64}, maintained_parameters::Array{Float64})
    index_of_cell_to_choose = rand(1: size(cluster,2)) # includes lower and upper number
    chosen_cell = cluster[:,index_of_cell_to_choose]
    random_number = rand(Float64)

    if (chosen_cell[parameter_dict["centre_1_position_x"]] == chosen_cell[parameter_dict["centre_2_position_x"]]) && (
                                    chosen_cell[parameter_dict["centre_1_position_y"]] == chosen_cell[parameter_dict["centre_2_position_y"]])
        #### I CELL PROTOCOL

        if random_number <= 200/2501
            chosen_cell = grow_i_cell(chosen_cell, maintained_parameters)
            #print("grow\n")
        else
            chosen_cell = migrate_i_cell(chosen_cell, maintained_parameters)
            #print("migrate")
        end
        #print("current=",chosen_cell[parameter_dict["current_radius"]], "max=",maintained_parameters[maintained_parameters_dict["R_max"]], "\n" )
        if chosen_cell[parameter_dict["current_radius"]] == maintained_parameters[maintained_parameters_dict["R_max"]]
            print("replace I with M")
            chosen_cell = replace_i_cell_with_m_cell(chosen_cell, maintained_parameters)
        end

    else
        #### M CELL PROTOCOL
        if random_number <=2506/2506
            chosen_cell=deform_m_cell(chosen_cell, maintained_parameters)
        elseif random_number <= 1/2506
            chosen_cell=rotate_m_cell(chosen_cell, maintained_parameters)
        else
            chosen_cell=migrate_m_cell(chosen_cell, maintained_parameters)
        end
        print("current def=",chosen_cell[parameter_dict["current_deformation"]] , "current radius=", chosen_cell[parameter_dict["current_radius"]], "\n" )
        #if chosen_cell[parameter_dict["current_deformation"]] >= (2*chosen_cell[parameter_dict["current_radius"]])+0.5
        if chosen_cell[parameter_dict["current_radius"]] <= maintained_parameters[maintained_parameters_dict["R_max"]]/(sqrt(2))

            print("Replace M with 2I \n")
            two_i_cells = replace_m_cell_with_two_i_cells(chosen_cell, maintained_parameters)
            #print("size of new cells=", size(two_i_cells), "\n")
            chosen_cell = two_i_cells[:,1]
            other_cell = two_i_cells[:,2]
            cluster = hcat(cluster, other_cell)
        end
    end
    cluster[:,index_of_cell_to_choose]=chosen_cell
    return cluster
end

function draw_an_i_cell(i_cell::Array{Float64})
    x_pos = i_cell[parameter_dict["centre_1_position_x"]]
    y_pos = i_cell[parameter_dict["centre_1_position_y"]]
    radius= i_cell[parameter_dict["current_radius"]]
    circle = pyplt.Circle((x_pos, y_pos), radius, alpha=0.1 , color=:green, lw=0.3)
    return circle
end
    #pyplt.plot(draw_a_circle(x_pos,y_pos,radius_for_draw), seriestype = [:shape,], lw=1, aspect_ratio = 1, linecolor = :black, fillalpha = 0)
function draw_an_m_cell(m_cell::Array{Float64})
    x_pos1 = m_cell[parameter_dict["centre_1_position_x"]]
    y_pos1 = m_cell[parameter_dict["centre_1_position_y"]]
    x_pos2 = m_cell[parameter_dict["centre_2_position_x"]]
    y_pos2 = m_cell[parameter_dict["centre_2_position_y"]]
    radius= m_cell[parameter_dict["current_radius"]]
    circle1 = pyplt.Circle((x_pos1, y_pos1), radius, alpha = 0.1 , color = :blue, lw=0.3)
    circle2 = pyplt.Circle((x_pos2, y_pos2), radius, alpha=0.1 , color = :blue, lw=0.3)
    circles::Array{} = [circle1, circle2]
    #pyplt.plot(draw_a_circle(x_pos1,y_pos1,radius_for_draw), seriestype = [:shape,], lw=1, aspect_ratio = 1, linecolor = :blue, fillalpha = 0)
    #pyplt.plot(draw_a_circle(x_pos2,y_pos2,radius_for_draw), seriestype = [:shape,], lw=1, aspect_ratio = 1, linecolor = :blue, fillalpha = 0)
    return circles
end


ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
using PyPlot
using PyCall
@pyimport matplotlib.pyplot as pyplt


function main()
    maintained_parameters = create_maintained_parameters(5.0, 0.05, 0.05, 0.000001, 0.05, 0.005)
    print(maintained_parameters[maintained_parameters_dict["R_max"]])
    cluster_parameters = create_cluster_parameters(5.0,1.0,1.0,2.0,2.0,5.0)
    starting_cell = create_i_cell(1.0,1.0,4.9, 0.0)
    cluster = initialise_cluster(starting_cell)
    print(cluster)
    i = 0
    while i <500
        #print("i=", i, "\n")
        cluster = make_a_move(cluster, maintained_parameters)
        #print(cluster, "\n")

        i += 1
    end
    print(transpose(cluster))
    print("there are", size(cluster,2), "cells\n")

    #plt = plot( bg = :white, xlim = (0.0, 20.0),  ylim = (-10.0, 10.0), framestyle = :none,
    #size = (400, 400),  legend = false,)
    #pyplt.figure(figsize = (5,5))
    pyplt.clf()
    ax=pyplt.gca()
    ax.set_xlim((-10, 10))
    ax.set_ylim((-10, 10))
    ax.set_aspect("equal")
    for i::Int16 in (1:1:size(cluster,2))
        print("i=", i, "\n")
        cell_to_draw = cluster[:,i]
        #radius_for_draw = cell_to_draw[parameter_dict["current_radius"]]
        if (cell_to_draw[parameter_dict["centre_1_position_x"]] == cell_to_draw[parameter_dict["centre_2_position_x"]]) && (
                                        cell_to_draw[parameter_dict["centre_1_position_y"]] == cell_to_draw[parameter_dict["centre_2_position_y"]])
            #I CELL PROTOCOL
            ax.add_artist(draw_an_i_cell(cell_to_draw))
        else
            circles = draw_an_m_cell(cell_to_draw)
            current_cell = circles[1]
            ax.add_artist(current_cell)
            current_cell = circles[2]
            ax.add_artist(current_cell)
            #end

        end
        fig = gcf()
        display(fig)
    end

end

main()
