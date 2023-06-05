% تابع محاسبهٔ f(x)
function result = f(x)
    result = x(1)^2 + 2 * x(2)^2 - 0.3 * cos(3 * pi * x(1)) - 0.4 * cos(4 * pi * x(2));
end

% حدس اولیه
x0 = [0.8; -0.25];

% پارامترهای الگوریتم
alpha = 0.075; % طول قدم
beta = 0.01; % عرض ناحیه عدم قطعیت

% متغیرهای مربوط به الگوریتم
x_prev = x0;
f_prev = f(x_prev);
iteration = 0;

% نمایش سرستون‌های جدول
disp('تعداد تکرار     x1       x2      f(x)');

% الگوریتم گرادیان
while true
    % محاسبهٔ گرادیان
    gradient = [2 * x_prev(1) + 0.3 * pi * sin(3 * pi * x_prev(1));
                4 * x_prev(2) + 0.4 * pi * sin(4 * pi * x_prev(2))];
    
    % محاسبهٔ جهت نزولی
    direction = -gradient;
    
    % جستجوی خط
    t = 1;
    while f(x_prev + t * direction) >= f_prev
        t = beta * t;
    end
    
    % به‌روزرسانی مقدار جدید
    x_new = x_prev + t * direction;
    f_new = f(x_new);
    
    % نمایش مقادیر در جدول
    disp(['iteration, x_prev', f_prev]);
    
    % شرایط خاتمه الگوریتم
    if norm(x_new - x_prev) < beta || abs(f_new - f_prev) < beta
        break;
    end
    
    % به‌روزرسانی مقادیر برای تکرار بعدی
    x_prev = x_new;
    f_prev = f_new;
    iteration = iteration + 1;
end

% نمایش مقدار بهینه و تابع در نقطه بهینه
disp('مقدار بهینهٔ x:');
disp(x_new);
disp('مقدار بهینهٔ f(x):');
disp(f_new);
