import numpy as np
from scipy.spatial.distance import euclidean as dist

class Coordinates:
    def __init__(self):
        self.raft_list = ['01','02','03',
                          '10','11','12','13','14',
                          '20','21','22','23','24',
                          '30','31','32','33','34',
                          '41','42','43']
        self.raftcentercoord = {}
        for raft in self.raft_list:
            tmpx = 127.*(int(raft[0]) - 2) # raft center x-coordinate in camera coordinate system
            tmpy = 127.*(int(raft[1]) - 2) # raft center y-coordinate in camera coordinate system
            self.raftcentercoord[raft] = [tmpx,tmpy]

        self.sensor_list = ['00','01','02',
                            '10','11','12',
                            '20','21','22']
        self.sensorcentercoord = {}
        for sensor in self.sensor_list:
            tmpx = 42.25*(int(sensor[0]) - 1) # sensor center x-coordinate in raft coordinate system
            tmpy = 42.25*(int(sensor[1]) - 1) # sensor center y-coordinate in raft coordinate system
            self.sensorcentercoord[sensor] = [tmpx,tmpy]
    
        
    def cam_to_raft(self,camcoord):
        '''
        For a given set of coordinates camcoord = (cam_x, cam_y) in the camera coordinate system, returns:
            - the ID of the raft whose center is the closest to camcoord
            - the corresponding coordinates in the raft frame, where raftcoord = (raft_x, raft_y) = (0,0) at the raft centre
        '''
        self.camcoord = camcoord
        self.dist2raft={}
        for raft in self.raft_list:
#            print(self.camcoord, self.raftcentercoord[raft])
            self.dist2raft[raft] = dist(camcoord, self.raftcentercoord[raft])
        self.raft_id = min(self.dist2raft, key=self.dist2raft.get)
        self.raftcoord = np.empty(2)
        self.raftcoord[0] = self.camcoord[0] - self.raftcentercoord[self.raft_id][0]
        self.raftcoord[1] = self.camcoord[1] - self.raftcentercoord[self.raft_id][1]
        
#        return self.raft_id, self.raftcoord
    
    def raft_to_sensor(self, raftcoord):
        '''
        For a given set of coordinates raftcoord = (raft_x, raft_y) in the raft coordinate system, returns:
            - the ID of the ccd whose center is the closest to raftcoord
            - the corresponding coordinates in the ccd frame, where (slotcoord_x, slotcoord_y) = (0,0) at the CCD centre
        '''
        self.raftcoord = raftcoord
        self.dist2sensor={}
        for sensor in self.sensor_list:
#            print(self.camcoord, self.raftcentercoord[raft])
            self.dist2sensor[sensor] = dist(raftcoord, self.sensorcentercoord[sensor])
        self.sensor_id = min(self.dist2sensor, key=self.dist2sensor.get)
        self.sensorcoord = np.empty(2)
        self.sensorcoord[0] = self.raftcoord[0] - self.sensorcentercoord[self.sensor_id][0]
        self.sensorcoord[1] = self.raftcoord[1] - self.sensorcentercoord[self.sensor_id][1]
        
        return self.sensor_id, self.sensorcoord

    def sensor_to_segment(self, sensorcoord):
        '''
        For a given set of coordinates (xcam, ycam) in the camera coordinate system, returns:
            - the ID of the raft whose center is the closest to (xcam, ycam)
            - the corresponding coordinates in the raft frame, where (slotcoord_x, slotcoord_y) = (0,0) at the CCD centre
        '''
        return self.seg_id, self.segcoord

########## And the other way around ############################# 
    
    def raft_to_cam(self, raft_id, raftcoord):
        '''
        For a given set of coordinates (cam_x, cam_y) in the camera coordinate system, returns:
            - the ID of the raft whose center is the closest to (xcam, ycam)
            - the corresponding coordinates in the raft frame, where (raft_x, raft_y) = (0,0) at the raft centre
        '''
        return self.camcoord
   
    def slot_to_raft(self, slot_id, slotcoord):
        '''
        For a given set of coordinates (raft_x, raft_y) in the camera coordinate system, returns:
            - the ID of the raft whose center is the closest to (xcam, ycam)
            - the corresponding coordinates in the raft frame, where (slotcoord_x, slotcoord_y) = (0,0) at the CCD centre
        '''
        return self.raft_id, self.raftcoord
    
    def segment_to_slot(self, seg_id, segcoord):
        '''
        For a given set of coordinates (xcam, ycam) in the camera coordinate system, returns:
            - the ID of the raft whose center is the closest to (xcam, ycam)
            - the corresponding coordinates in the raft frame, where (slotcoord_x, slotcoord_y) = (0,0) at the CCD centre
        '''
        return self.slot_id, self.slotcoord
    
    
