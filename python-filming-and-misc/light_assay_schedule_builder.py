import random
import datetime

# User inputs the first intended filming day, where fish is 17 dpf
print("Enter the date at which the brood is at 17 dpf");
yy = int(input("Year: "));
mm = int(input("Month (int): "));
dd = int(input("Day: "));

# Create datetime object
date = datetime.date(yy, mm, dd);
   
# value of 0 indicates 17 dpf
# Each value of these lists is a day-number offset
start_days = [0, 4, 8, 12];
end_days = [1, 5, 9, 13];
    
# Each group has 5 fish from the same brood, and represents one lighting condition
group1list = [];
group2list = [];
    
# Each treatment can only be the starting condition on 2 sets of consecutive days
for i in range(2):
    
    # Randomly assigns starting days
    i = random.randint(0, len(start_days) - 1);
    group1list.append(start_days[i]);
    group2list.append(end_days[i]);
    
    # Avoids repeat dates
    start_days.remove(start_days[i]);
    end_days.remove(end_days[i]);

# Fills in the rest of the days
while len(start_days) > 0:
    group1list.append(end_days[0]);
    group2list.append(start_days[0]);
    start_days.remove(start_days[0]);
    end_days.remove(end_days[0]);
    
group1list.sort();
group2list.sort();

#print(group1list);
#print(group2list);

print("\nLight treatment condition is recorded on the following days:");
for i in group1list:
    print(f'Record on {date + datetime.timedelta(i)} at {17+i} dpf');
    #print(f'At {17+i} dpf\n');

print("\nDark treatment condition is recorded on the following days:");
for i in group2list:
    print(f'Record on {date + datetime.timedelta(i)} at {17+i} dpf');
    #print(f'At {17+i} dpf\n');
