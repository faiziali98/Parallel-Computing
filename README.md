# Cardiac Electrophysiology Simulation

I have added code for 1D and 2D version of the assignment in files:

1. cardiacsimSerial.C (Serial)
2. cardiacsim1D.C (1D Version)
3. cardiacsim2D.C (2D Version)

You can run the code using commands:

```bash
make
mpirun -n (number of threads) ./cardiacsim2D -n (grid size) -t (time step)
mpirun -n (number of threads) ./cardiacsim1D -n (grid size) -t (time step)
./cardiacsimSerial -n (grid size) -t (time step)
```

Last but not least, we still have to:

1. Use openMP
2. Optamize the code according to the last few lectures
3. Test the code
