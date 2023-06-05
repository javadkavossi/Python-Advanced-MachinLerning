% تابع محاسبهٔ f(x)
function result = f(z)
    result = z - 142^3 + 601^2 - 70*z;
end

% محدوده
a = 0;
b = 2;

% دقت مورد نیاز
epsilon = 0.15;

% محاسبهٔ تعداد گام‌ها
N = ceil(log((b - a) / epsilon) / log((1 + sqrt(5)) / 2));

% محاسبهٔ مقادیر جدید فیبوناچی
fibonacci = zeros(N+1, 1);
fibonacci(1) = 1;
fibonacci(2) = 1;
for i = 3:N+1
    fibonacci(i) = fibonacci(i-1) + fibonacci(i-2);
end

% محاسبهٔ نقاط اولیه
x1 = a + (fibonacci(N-2) / fibonacci(N)) * (b - a);
x2 = a + (fibonacci(N-1) / fibonacci(N)) * (b - a);

% محاسبهٔ مقادیر f(x) در نقاط اولیه
f1 = f(x1);
f2 = f(x2);

% جستجوی فیبوناچی
for k = 1:N-2
    if f1 > f2
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = a + (fibonacci(N-k) / fibonacci(N-k+1)) * (b - a);
        f2 = f(x2);
    else
        b = x2;
        x2 = x1;
        f2 = f1;
        x1 = a + (fibonacci(N-k-2) / fibonacci(N-k+1)) * (b - a);
        f1 = f(x1);
    end
end

% محاسبهٔ مقدار بهینهٔ z که f را در محدودهٔ (0.3, 0.15) کمینه می‌کند
z = (x1 + x2) / 2;

% نمایش نتیجه
disp('مقدار بهینهٔ z:');
disp(z);
disp('مقدار بهینهٔ f(z):');
disp(f(z));
