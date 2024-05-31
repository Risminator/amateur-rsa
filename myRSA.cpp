#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <bitset>
#include <string>

// Проверка числа на простоту тестом Миллера-Рабина
bool MillerRabinTest(uint64_t);
// Генерация 32-битного числа
uint64_t generateNumber();
// Генерация пары двух простых 32-битных чисел
std::vector<uint64_t> generatePrimePair();
// НОД
uint64_t gcd(uint64_t, uint64_t);
// Мультипликативное обратное
uint64_t mult_inverse(uint64_t, uint64_t);
// Умножение по модулю
uint64_t mulMod(uint64_t, uint64_t, uint64_t);
// Быстрое возведение в степень по модулю
uint64_t powMod(uint64_t, uint64_t, uint64_t);

int main() {
    std::cout << "Простые числа p и q:\t";
	std::vector<uint64_t> PrimeNumbers = generatePrimePair();
	for (int i = 0; i < PrimeNumbers.size(); i++) {
		std::cout << PrimeNumbers[i] << ' ';
	}
	std::cout << "\n";
    
	uint64_t N = PrimeNumbers[0] * PrimeNumbers[1];
	std::cout << "Модуль N = p * q:\t" << N << "\n";

	uint64_t euler = (PrimeNumbers[0] - 1) * (PrimeNumbers[1] - 1);
	std::cout << "Функция Эйлера phi(N):\t" << euler << "\n";

	// Генерируем достаточно большое e для защиты от атаки Хастада
	int e;
	do {
		e = 100000 + rand() % 900000;
	} while (gcd(e, euler) != 1 || e >= euler);	
	std::cout << "Открытая экспонента:\t" << e << "\n";

	uint64_t d = mult_inverse(e, euler);
	std::cout << "Секретная экспонента d:\t" << d << "\n";

	std::cout << "Открытый ключ (e, N):\t(" << e << ", " << N << ")\n";
	std::cout << "Закрытый ключ (d, N):\t(" << d << ", " << N << ")\n";

	std::cout << "Введите сообщение\n";
	std::string message = "";
	std::getline(std::cin, message);

	uint64_t* encryptedMessage = new uint64_t[message.size() + 1];

	// При шифровании используется бинарное возведение в степень
	std::cout << "Зашифрованное сообщение: ";
	for (int i = 0; i < message.size(); i++) {
		encryptedMessage[i] = powMod((uint64_t)message[i], e, N);
		std::cout << encryptedMessage[i] << " ";
	}
	std::cout << "\n";

	// Дешифруем сообщение с использованием китайской теоремы об остатках
	// Также используется бинарное возведение в степень
	std::cout << "Дешифрованное сообщение: ";
	uint64_t dP = d % (PrimeNumbers[0] - 1);
	uint64_t dQ = d % (PrimeNumbers[1] - 1);
	uint64_t qInv = mult_inverse(PrimeNumbers[1], PrimeNumbers[0]);
	for (int i = 0; i < message.size(); i++) {
		uint64_t m1 = powMod(encryptedMessage[i], dP, PrimeNumbers[0]);
		uint64_t m2 = powMod(encryptedMessage[i], dQ, PrimeNumbers[1]);
		uint64_t h = qInv * (m1 - m2) % PrimeNumbers[0];
		uint64_t m = (m2 + h * PrimeNumbers[1]) % N;
		std::cout << (char)m;
	}
	std::cout << "\n";
	delete[] encryptedMessage;
}

// Функция умножения числа a на b по модулю n
uint64_t mulMod(uint64_t a, uint64_t b, uint64_t n) {
	uint64_t res = 0;

	while (a != 0) {
		if (a & 1) {
			res = (res + b) % n;
		}
		a >>= 1;
		b = (b << 1) % n;
	}
	return res;
}

// Функция возведения числа a в степень b по модулю n
uint64_t powMod(uint64_t a, uint64_t b, uint64_t n) {
	uint64_t x = 1;

	a %= n;

	while (b > 0) {
		if (b & 1) {
			x = mulMod(x, a, n); 
		}
		a = mulMod(a, a, n);
		b >>= 1;
	}
	return x % n;
}

// Функция для проверки числа на простоту с помощью теста Миллера-Рабина
bool MillerRabinTest(uint64_t n) {
    int k = 20;
    // Обработка исключений: проверка на простоту тривиальных случаев
    if (n <= 1)
        return false;
    if (n <= 3)
        return true;

    // Находим целое число t и нечетное число u, такие что n-1 = (2^t)*u
    long long u = n - 1;
    int t = 0;
    while (u % 2 == 0) {
        u = u / 2;
        t++;
    }

    // Повторяем тест k раз
    for (int i = 0; i < k; i++) {
        // Выбираем случайное целое число a в интервале [2, n-2]
        uint64_t a = 2 + rand() % (n - 3);

        // Вычисляем a^u % n
        uint64_t x = powMod(a, u, n);

        // Если результат x равен 1 или n-1, то число вероятно простое, проверяем следующее
        if (x == 1 || x == n - 1)
            continue;

        // Применяем алгоритм Миллера-Рабина
        for (int j = 0; j < t - 1; j++) {
            x = powMod(x, 2, n);
            // Если x равно 1, число составное
            if (x == 1)
                return false;
            // Если x равно n-1, число вероятно простое, проверяем следующее
            if (x == n - 1)
                break;
        }
        // Если не выполнились условия внутреннего цикла, число составное
        if (x != n - 1)
            return false;
    }
    // Если не выполнились условия внешнего цикла, число вероятно простое
    return true;
}

// Функция для генерации большого  числа 
uint64_t generateNumber() {
    constexpr uint32_t bits = 31;
    std::bitset<bits> a;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<short> distr(0, 1);

    for (uint32_t i = 0; i < bits; i++) {
        a[i] = distr(gen);
    }

    a[0] = 1;
    a[bits - 1] = 1;

    return a.to_ullong();
}

// Функция для генерации двух простых чисел
std::vector<uint64_t> generatePrimePair() {
    uint64_t num1, num2;
    do {
        num1 = generateNumber();
        // Sleep(2);
        num2 = generateNumber();
    } while ((MillerRabinTest(num1) && MillerRabinTest(num2) && num1 != num2) != 1);
    std::vector<uint64_t> PrimeNumbers = { num1, num2 };
    return PrimeNumbers;
}

// Функция для нахождения НОД двух чисел
uint64_t gcd(uint64_t a, uint64_t b) {
	// Если одно из чисел равно 0, то НОД равен другому числу
	if (a == 0)
		return b;
	if (b == 0)
		return a;

	// Применяем алгоритм Евклида для нахождения НОД
	while (b != 0) {
		uint64_t remainder = a % b;
		a = b;
		b = remainder;
	}
	return a;
}

std::tuple<uint64_t, uint64_t> sequence(uint64_t a, uint64_t b, uint64_t s0, uint64_t s1) {
	uint64_t counter = 0;
	while (b != 0) {
		++counter;
		uint64_t q = a / b;
		uint64_t tmp = b;
		b = a % b;
		a = tmp;
		tmp = s1;
		s1 = q * s1 + s0;
		s0 = tmp;
	}
	return std::make_tuple(s0, counter);
}

uint64_t mult_inverse(uint64_t a, uint64_t n) {
	a = a % n;
	if (gcd(a, n) != 1) {
		std::cout << "Error in mult_inverse with a = " << a << ", n = " << n << '\n';
		return 0;
	}
	std::tuple<uint64_t, uint64_t> xsize = sequence(a, n, 1, 0);

	int coef;
	if (!(std::get<1>(xsize) % 2)) coef = 1;
	else coef = -1;

	return (coef == 1 ? std::get<0>(xsize) : n - std::get<0>(xsize));
}