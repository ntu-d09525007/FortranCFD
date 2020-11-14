from paraview.simple import *

timesteps = GetTimeKeeper().TimestepValues
view = GetAnimationScene()
renderView1 = GetActiveViewOrCreate('RenderView')

for i in range(len(timesteps)):
	view.AnimationTime = timesteps[i]
	ExportView('D:/RG30/' +str(i)+'.pov', view=renderView1)
