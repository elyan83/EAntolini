import subprocess

subprocess.call(["echo", "I am doing something."])

print("Sleeping for ten seconds using subprocess.call.")
subprocess.call(["sleep", "10"])
print("I have just run the sleep command.")
subprocess.call(["echo", "I am doing something else. I had to wait for sleep to complete."])

print("Sleeping for ten seconds using subprocess.Popen.")
sleepproc=subprocess.Popen(["sleep", "10"])
print("I have just called the sleep command using Popen.")
print("I'm going to run something while sleep is running.")
subprocess.call(["echo", "I am doing something else while sleep is running."])
print("Now I'm going to wait for sleep to finish.")
sleepproc.wait()
print("All done.")

