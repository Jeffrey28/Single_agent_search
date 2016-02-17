#!/home/mpc/.virtualenvs/research/bin/python

import rospy
import json
import requests
import time
from std_msgs.msg import *
from geometry_msgs.msg import PoseStamped, Vector3, Pose2D


topic = 'mpc_solver_feed'

condition = Float64MultiArray()

def talker():
    pub = rospy.Publisher(topic, Float64MultiArray, queue_size=10)
    rospy.init_node('condition_provider')
    r = rospy.Rate(5)
    condition.data = [0.0, 3.3, 4.4, 5.5]
    pub.publish(condition)
    time.sleep(2)


if __name__ == '__main__':
    try:
        talker()
    except rospy.ROSInterruptException: pass
