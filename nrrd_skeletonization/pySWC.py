# -*- coding: utf-8 -*-
"""
pyswc is designed for the manipulation of swc files
"""
import numpy as np

class Swc:
    
    #indices for elements in the swc file
    node_number = 0
    node_label = 1
    node_x = 2
    node_y = 3
    node_z = 4
    node_radius = 5
    node_parent = 6
    
    def __init__(self, filename):
        
        self.data = np.genfromtxt(filename)
        self.num_nodes = self.data.shape[0]
        self.surface_area = self.find_surface_area()
        
    def save_file(self, filename):
        
        np.savetxt(filename, self.data, fmt='%i %i %4.12f %4.12f %4.12f %4.5f %i')  

    def relabel(self, old, new):
        
        label = self.node_label
        
        for node in range(0,self.num_nodes):
            if self.data[node, label] == old:
                self.data[node, label] = new

    def scale(self, scale_factor, scale_radii=True):
        
        x, y, z, radius = self.node_x, self.node_y, self.node_z, self.node_radius
        
        for node in range(0,self.num_nodes):
            
            self.data[node, x] = self.data[node, x] * scale_factor
            self.data[node, y] = self.data[node, y] * scale_factor
            self.data[node, z] = self.data[node, z] * scale_factor
            if scale_radii:
                self.data[node, radius] = self.data[node, radius] * scale_factor
            
        self.surface_area =  self.find_surface_area()
        
    def translate(self, x_translation, y_translation, z_translation):
        
        x, y, z = self.node_x, self.node_y, self.node_z
        
        for node in range(0, self.num_nodes):
            
            self.data[node, x] = self.data[node, x] + x_translation
            self.data[node, y] = self.data[node, y] + y_translation
            self.data[node, z] = self.data[node, z] + z_translation
            
    def find_surface_area(self):
        """
        Models the swc file as a series of conical frustums to calculate Surface area
        """
        x, y, z, radius = self.node_x, self.node_y, self.node_z, self.node_radius
        
        total_area = 0.0
        for node in range(0, self.num_nodes-1):
            x_dist = np.abs(self.data[node, x] - self.data[node+1, x])
            y_dist = np.abs(self.data[node, y] - self.data[node+1, y])
            z_dist = np.abs(self.data[node, z] - self.data[node+1, z])
            total_dist = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)
            section_area = (np.pi*(self.data[node, radius]+self.data[node+1, radius]))*(np.sqrt(total_dist**2+(self.data[node+1, radius]-self.data[node, radius])**2))
            total_area += section_area
        
        return(total_area)

    def adjust_surface_area(self, goal_surface_area, error_rate=0.001):
        """
        adjusts the surface area of the swc file to match a goal value.
        This is done iteratively by scaling the radii, checking the surface area
        and then rescaling the radii until the surface area matches the goal area
        within the error rate (default=0.001)
        """
        radius = self.node_radius
        
        while np.abs(self.surface_area - goal_surface_area)/goal_surface_area > error_rate:
            new_swc = self.data
            scaling_factor = goal_surface_area / self.surface_area
            for node in range(0, self.num_nodes):
                new_swc[node, radius] = scaling_factor * self.data[node, radius]
            self.data = new_swc
            self.surface_area = self.find_surface_area()
    
    def rotate(self, x_rotation_deg, y_rotation_deg, z_rotation_deg):
        import math
        
        x, y, z = self.node_x, self.node_y, self.node_z

        x_rotation = x_rotation_deg * math.pi / 180
        y_rotation = y_rotation_deg * math.pi / 180
        z_rotation = z_rotation_deg * math.pi / 180
        
        for node in range(0,self.num_nodes):

            #x rotations
            self.data[node,x] = self.data[node,x]
            self.data[node,y] = self.data[node,y] * math.cos(x_rotation) - self.data[node,z] * math.sin(x_rotation)
            self.data[node,z] = self.data[node,y] * math.sin(x_rotation) + self.data[node,z] * math.cos(x_rotation) 
    
            #y rotations
            self.data[node,x] = self.data[node,x] * math.cos(y_rotation) + self.data[node,z] * math.sin(y_rotation)
            self.data[node,y] = self.data[node,y]
            self.data[node,z] = self.data[node,z] * math.cos(y_rotation) - self.data[node,x] * math.sin(y_rotation)
        
            #z rotations
            self.data[node,x] = self.data[node,x] * math.cos(z_rotation) - self.data[node,y] * math.sin(z_rotation)
            self.data[node,y] = self.data[node,x] * math.sin(z_rotation) + self.data[node,y] * math.cos(z_rotation)
            self.data[node,z] = self.data[node,z]

    def find_children(self, node):

        number, parent = self.node_number, self.node_parent
        children = []

        for other_node in range(0,self.num_nodes):

            if self.data[other_node, parent] == self.data[node, number]:
                children.append(other_node)

        return(children)

    def halve_num_nodes(self):

        number, parent = self.node_number, self.node_parent
        delete = 0
        for node in range(0, self.num_nodes):
            children = self.find_children(node)
            if len(children) == 1 and delete == 0:

                current_parent = self.data[node, parent]
                self.data[children[0], parent] = current_parent
                del self.data[node, :]

                for rem_node in range(node, self.num_nodes):
                    self.data[rem_node, number] -= 1
                    self.data[rem_node, parent] -= 1

                delete += 1

            delete = delete % 2
